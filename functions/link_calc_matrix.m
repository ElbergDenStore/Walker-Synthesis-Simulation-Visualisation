function Link_2D = link_calc_matrix(el_mat, range_mat, lat_vec, lon_vec, link_cfg, general_config)
    % Inputs are now NumUEs x nT matrices
    Link_2D.Frequency = link_cfg.f;
    Link_2D.Bandwidth = link_cfg.B;
    
    % 1. FSPL (Matrix Math)
    c = physconst('LightSpeed');
    lambda = c / link_cfg.f;
    Link_2D.FSPL = 10 * log10(((4 * pi * range_mat) / lambda).^2);
    
    % 2. Steering Loss (Matrix Math)
    [Link_2D.Rx_steering_loss, Link_2D.Tx_steering_loss] = Steering_loss_calc_2D(el_mat, range_mat, link_cfg.Tx_type, link_cfg.Rx_type);
    
    % 3. Absorption (Matrix Math for Fast Mode)
    [Link_2D.Absorption_At, Link_2D.T_antenna] = Absorption_calc_2D(link_cfg.f, el_mat, lat_vec, lon_vec, general_config);
    
    % 4. Total Loss
    Link_2D.Total_loss = Link_2D.FSPL + Link_2D.Absorption_At + Link_2D.Rx_steering_loss + Link_2D.Tx_steering_loss;
    
    % 5. Adjusted Power & PFD
    if (isfield(link_cfg, 'Direction') && link_cfg.Direction == "UL")
        Link_2D.T_antenna = zeros(size(el_mat)) + 290;
        Link_2D.Adjusted_EIRP_dBm = zeros(size(el_mat)) + link_cfg.Max_EIRP_dBm;
        Link_2D.PFD_W_MHz = nan(size(el_mat)); 
    else
        Link_2D.Adjusted_EIRP_dBm = Adjust_tx_power_2D(el_mat, range_mat, general_config.Min_elevation_UE, general_config.Orbit_height, link_cfg.Tx_type, link_cfg.Max_EIRP_dBm); 
        
        area = 10*log10(4*pi.*range_mat.^2); 
        Link_2D.PFD_W_MHz = Link_2D.Adjusted_EIRP_dBm - 30 - area - Link_2D.Tx_steering_loss - 10*log10(link_cfg.B/1e6);
    end
    
    % 6. SNR & Throughput (Matrix Math)
    F_lin = 10^(link_cfg.NF/10);
    T_rx = (F_lin - 1) * 290;
    T_sys = Link_2D.T_antenna + T_rx;
    Link_2D.P_noise = (10*log10(T_sys) + 10*log10(1.38e-23)) + 30 + 10*log10(link_cfg.B);
    
    Link_2D.Rx_Power = Link_2D.Adjusted_EIRP_dBm + link_cfg.G_rx - Link_2D.Total_loss;
    Link_2D.SNR = Link_2D.Rx_Power - Link_2D.P_noise;
    
    % Throughput calculation
    Link_2D.Throughput = zeros(size(Link_2D.SNR));
    valid_idx = ~isnan(Link_2D.SNR);
    
    SNR_lin = 10.^(Link_2D.SNR(valid_idx)/10);
    SNR_eff = 10^(2/10); % 2dB effective SNR assumed from your original code
    Link_2D.Throughput(valid_idx) = link_cfg.B * 0.56 * 1 .* log2(1 + SNR_lin ./ SNR_eff);
end

% --- Subfunctions adjusted for matrices ---
function [loss_rx, loss_tx] = Steering_loss_calc_2D(el_mat, range_mat, tx_type, rx_type)
    Re = 6378.14e3;     
    theta_rx = 90 - el_mat; 
    R_sat = sqrt(Re^2 + range_mat.^2 + 2 * Re .* range_mat .* sind(el_mat));
    theta_tx = asind((Re .* cosd(el_mat)) ./ R_sat); 
    
    loss_tx = zeros(size(el_mat)); 
    loss_rx = zeros(size(el_mat));
    if contains(tx_type, 'array'), loss_tx = -10 * log10(cosd(theta_tx).^1.5); end
    if contains(rx_type, 'array'), loss_rx = -10 * log10(cosd(theta_rx).^1.5); end
end

function adjusted_tx_power = Adjust_tx_power_2D(el_mat, range_mat, min_elev, orb_ht, tx_type, max_pwr)
    Re = 6378.14e3;     
    worst_sin = (Re * cosd(min_elev)) / (orb_ht + Re);
    worst_dist_inc = (-Re*sind(min_elev)+sqrt(Re^2*sind(min_elev)^2 -(Re^2-(Re+orb_ht)^2)))/orb_ht; 
    
    dist_dec = -20*log10(range_mat ./ (worst_dist_inc*orb_ht)); 
    
    if contains(tx_type, 'array')
        worst_beam_inc = 1/(cosd(asind(worst_sin)).^1.5);
        theta_tx = asind((Re .* cosd(el_mat)) ./ (orb_ht + Re));
        beam_dec = -10*log10((1./(cosd(theta_tx).^1.5))/worst_beam_inc);
        adjusted_tx_power = max_pwr - dist_dec - beam_dec;
    else
        adjusted_tx_power = max_pwr - dist_dec;
    end
end

function [At_mat, Tsky_mat] = Absorption_calc_2D(f, el_mat, lat_vec, lon_vec, Cfg)
    % Fast Mode: Instantaneous matrix fill
    if isfield(Cfg, 'Use_P618') && Cfg.Use_P618 == false
        static_loss = 0.5;
        if isfield(Cfg, 'Simple_Atmospheric_Loss_dB'), static_loss = Cfg.Simple_Atmospheric_Loss_dB; end
        At_mat = zeros(size(el_mat)) + static_loss;
        Tsky_mat = zeros(size(el_mat)) + 290;
        return;
    end
    
    % Precision Mode (P.618)
    % Because P618 relies on LAT/LON specific lookups, we MUST loop over UEs here.
    % But we only loop over the spatial dimension, keeping time vectorized.
    At_mat = zeros(size(el_mat));
    Tsky_mat = zeros(size(el_mat));
    el_grid = 20:10:90; 
    
    num_pts = length(lat_vec);
    
    use_par = false;
    if isfield(Cfg, 'use_parallel') && Cfg.use_parallel
        use_par = true;
    end
    
    if ~isempty(getCurrentTask())
        use_par = false; % Cannot use DataQueue/parfor on a worker
    end
    
    if use_par
        updateLiveScriptProgress(num_pts, true);
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(~) updateLiveScriptProgress(num_pts, false));
        
        parfor u = 1:num_pts
            grid_At = zeros(1, length(el_grid)); grid_Tsky = zeros(1, length(el_grid));
            ws = warning('off', 'all'); 
            for i = 1:length(el_grid)
                link_cfg = p618Config('Frequency',f,'ElevationAngle',el_grid(i),'Latitude',lat_vec(u),'Longitude',lon_vec(u),'TotalAnnualExceedance',1);      
                [pl, ~, tsky] = p618PropagationLosses(link_cfg, 'StationHeight', 0);
                grid_At(i) = pl.At; grid_Tsky(i) = tsky;
            end
            warning(ws); 
            
            At_mat(u, :) = interp1(el_grid, grid_At, el_mat(u,:), 'linear');
            Tsky_mat(u, :) = interp1(el_grid, grid_Tsky, el_mat(u,:), 'linear');
            
            send(dq, []);
        end
    else
        old_warn = warning('off', 'all'); 
        for u = 1:num_pts
            grid_At = zeros(1, length(el_grid)); grid_Tsky = zeros(1, length(el_grid));
            for i = 1:length(el_grid)
                link_cfg = p618Config('Frequency',f,'ElevationAngle',el_grid(i),'Latitude',lat_vec(u),'Longitude',lon_vec(u),'TotalAnnualExceedance',1);      
                [pl, ~, tsky] = p618PropagationLosses(link_cfg, 'StationHeight', 0);
                grid_At(i) = pl.At; grid_Tsky(i) = tsky;
            end
            At_mat(u, :) = interp1(el_grid, grid_At, el_mat(u,:), 'linear');
            Tsky_mat(u, :) = interp1(el_grid, grid_Tsky, el_mat(u,:), 'linear');
            
            if ~use_par && isempty(getCurrentTask())
                updateLiveScriptProgress(num_pts, false);
            end
        end
        warning(old_warn); 
    end
end