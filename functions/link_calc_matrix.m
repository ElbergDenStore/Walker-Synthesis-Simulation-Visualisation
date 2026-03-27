function Link_2D = link_calc_matrix(el_mat, az_mat, range_mat, lat_vec, lon_vec, link_cfg, general_config)
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

    if (isfield(link_cfg, 'Direction') && link_cfg.Direction == "UL")
        BeamGrid = calculate_Beams(link_cfg.f, link_cfg.G_rx, general_config.Orbit_height, general_config.Min_elevation_UE, general_config.FRF);
    else
        BeamGrid = calculate_Beams(link_cfg.f, link_cfg.G_tx, general_config.Orbit_height, general_config.Min_elevation_UE, general_config.FRF);
    end
    
    [sir_lin_mat, mb_idx_mat] = Interference_calc_2D(el_mat, az_mat, BeamGrid, general_config);
    Link_2D.SIR = 10 * log10(sir_lin_mat);
    
    % Convert signal and noise back to linear milliwatts
    S_mW = 10.^(Link_2D.Rx_Power / 10);
    N_mW = 10.^(Link_2D.P_noise / 10);
    SIR_lin = 10.^(Link_2D.SIR / 10);
    
    % Calculate exact Interference power
    I_mW = S_mW ./ SIR_lin;
    
    % Final SINR
    SINR_lin = S_mW ./ (I_mW + N_mW);
    Link_2D.SINR = 10 * log10(SINR_lin);

    % Throughput calculation
    Link_2D.Throughput = zeros(size(Link_2D.SNR));
    valid_idx = ~isnan(Link_2D.SNR);
    
    SINR_lin_valid = 10.^(Link_2D.SINR(valid_idx)/10);
    %%%% Modified shannon %%%
    BW_eff = 0.56;
    eta = 1;
    SNR_eff_dB = 2;
    SNR_eff = 10^(SNR_eff_dB/10);

    % Calculate the raw capacity of the beam
    raw_throughput = link_cfg.B * BW_eff * eta .* log2(1 + SINR_lin_valid ./ SNR_eff);
    %%%% modified shannon end %%%
    
    % Check if the Bandwidth Sharing flag is turned on
    if isfield(general_config, 'Share_bandwidth') && general_config.Share_bandwidth
        NumUEs = size(mb_idx_mat, 1);
        nT = size(mb_idx_mat, 2);
        
        % Create a time index matrix the exact same size as the UEs
        time_matrix = repmat(1:nT, NumUEs, 1);
        
        valid_mask = ~isnan(mb_idx_mat);
        mb_valid = mb_idx_mat(valid_mask);
        time_valid = time_matrix(valid_mask);
        
        % Create a unique ID for every single "Beam at Time T" combination
        unique_beam_time_ids = mb_valid + (time_valid * BeamGrid.num_beams);
        
        % Vectorized counting: How many UEs share the exact same Beam-Time ID?
        [~, ~, ic] = unique(unique_beam_time_ids);
        counts = accumarray(ic, 1);
        users_sharing_beam = counts(ic);
        
        % Re-inflate the counts back into a 2D matrix
        users_per_beam_mat = ones(NumUEs, nT);
        users_per_beam_mat(valid_mask) = users_sharing_beam;
        
        % Divide the throughput by the number of users sharing the beam!
        Link_2D.Throughput(valid_idx) = raw_throughput ./ users_per_beam_mat(valid_idx);
    else
        % No sharing, every user gets 100% of the beam's capacity
        Link_2D.Throughput(valid_idx) = raw_throughput;
    end
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

    dq = parallel.pool.DataQueue;
    updateLiveScriptProgress(Cfg.NumUEs, true); 
    afterEach(dq, @(~) updateLiveScriptProgress(Cfg.NumUEs, false));
    
    old_warn = warning('off', 'all'); 
    parfor (idx = 1:Cfg.NumUEs, Cfg.Num_workers)
        grid_At = zeros(1, length(el_grid)); grid_Tsky = zeros(1, length(el_grid));
        for i = 1:length(el_grid)
            link_cfg = p618Config('Frequency',f,'ElevationAngle',el_grid(i),'Latitude',lat_vec(idx),'Longitude',lon_vec(idx),'TotalAnnualExceedance',1);      
            [pl, ~, tsky] = p618PropagationLosses(link_cfg, 'StationHeight', 0);
            grid_At(i) = pl.At; grid_Tsky(i) = tsky;
        end
        % Interpolate for this specific UE's time series
        At_mat(idx, :) = interp1(el_grid, grid_At, el_mat(idx,:), 'linear');
        Tsky_mat(idx, :) = interp1(el_grid, grid_Tsky, el_mat(idx,:), 'linear');
        send(dq, []);
    end
    warning(old_warn); 
end