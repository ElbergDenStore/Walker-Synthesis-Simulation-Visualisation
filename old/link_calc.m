function Link = link_calc(UE, link_cfg, general_config)
    el_vec    = [UE.SimData.Elevation_deg];
    range_vec = [UE.SimData.Range];
    lat       = UE.Lat;
    lon       = UE.Lon;
    Link.Frequency = link_cfg.f;
    Link.Bandwidth = link_cfg.B;
    
    Link.FSPL = FSPL_calc(link_cfg.f, range_vec);
    [Link.Absorption, Link.T_antenna] = Absorption_calc(link_cfg.f, el_vec, lat, lon);
    [Link.Rx_steering_loss, Link.Tx_steering_loss] = Steering_loss_calc(el_vec, range_vec, link_cfg.Tx_type, link_cfg.Rx_type);
    Link.Total_loss = Link.FSPL + Link.Absorption.At + Link.Rx_steering_loss + Link.Tx_steering_loss;
    
    if (isfield(link_cfg, 'Direction') && link_cfg.Direction == "UL")
        Link.T_antenna = zeros(size(Link.T_antenna)) + 290;
        Link.Adjusted_EIRP_dBm = link_cfg.Max_EIRP_dBm;
    else
        Link.Adjusted_EIRP_dBm = Adjust_tx_power(el_vec, range_vec, general_config.Min_elevation_UE, general_config.Orbit_height, link_cfg.Tx_type, link_cfg.Max_EIRP_dBm); % only for DL
        Link.PFD_W_MHz = pfd_calc(Link.Adjusted_EIRP_dBm, Link.Bandwidth, Link.Tx_steering_loss, range_vec);
    end

    Link.P_noise = Noise_density_calc(Link.T_antenna, link_cfg.NF) + 10*log10(link_cfg.B);
    Link.Rx_Power = Link.Adjusted_EIRP_dBm + link_cfg.G_rx - Link.Total_loss;
    Link.SNR = Link.Rx_Power - Link.P_noise;
    
    Link.Throughput = nan(size(Link.SNR));
    valid_idx = ~isnan(Link.SNR);
    Link.Throughput(valid_idx) = Throughput_calc(Link.SNR(valid_idx), link_cfg.B);
    Link.Throughput(~valid_idx) = 0; % Force 0 for NaN
end

function pfd_dBW_MHz = pfd_calc(adjusted_EIRP_dBm, bandwidth, steer_loss_tx, range_vec)
    area = 10*log10(4*pi.*range_vec.^2); % area of sphere
    pfd_dBW_MHz = adjusted_EIRP_dBm - 30 - area - steer_loss_tx - 10*log10(bandwidth/1e6); %subtract 30 to get in Watts
end

function adjusted_tx_power = Adjust_tx_power(el_vec, range_vec, min_elevation, orbit_height, tx_type, max_power)
    Re = 6378.14e3;     
    cos_exponent = 1.5; 
    
    % Worst case
    worst_sin_theta_tx = (Re * cosd(min_elevation)) / (orbit_height + Re);
    worst_theta_tx =  asind(worst_sin_theta_tx);
    worst_beamHPBW_increase = 1/(cosd(worst_theta_tx).^cos_exponent);
    worst_distance_increase = (-Re*sind(min_elevation)+sqrt(Re^2*sind(min_elevation)^2 -(Re^2-(Re+orbit_height)^2)))/orbit_height; % from law of cosines + trust me bro equation manipulation


    sin_theta_tx = (Re .* cosd(el_vec)) ./ (orbit_height + Re);
    theta_tx = asind(sin_theta_tx);
    
    distance_decrease = -20*log10(range_vec ./ (worst_distance_increase*orbit_height)); % 20*log10 as power decreases by range^2

    if contains(tx_type, 'array')
        beamHPBW_increase = 1./(cosd(theta_tx).^cos_exponent);
        beamHPBW_decrease = -10*log10(beamHPBW_increase/worst_beamHPBW_increase);
        adjusted_tx_power = max_power - distance_decrease - beamHPBW_decrease;
    else
        adjusted_tx_power = max_power - distance_decrease;
    end
end

function loss_dB = FSPL_calc(f, range_vec)
    c = physconst('LightSpeed');
    lambda = c / f;
    loss_dB = 10 * log10(((4 * pi * range_vec) / lambda).^2);
end

function Throughput = Throughput_calc(SNR_dB,B)
    BW_eff = 0.56;
    eta = 1;
    SNR_eff_dB = 2;
    SNR_eff = 10^(SNR_eff_dB/10);
    SNR = 10.^(SNR_dB/10);
    Throughput = B*BW_eff*eta.*log2(1+SNR./SNR_eff);
end

function [loss_struct, sky_temp_K] = Absorption_calc(f, el_vec, lat, lon)
    el_grid = 20:10:90; 
    n_grid = length(el_grid);
    grid_Ag = zeros(1, n_grid); grid_Ac = zeros(1, n_grid);
    grid_Ar = zeros(1, n_grid); grid_As = zeros(1, n_grid);
    grid_At = zeros(1, n_grid); grid_Tsky = zeros(1, n_grid);
    
    old_warn_state = warning('off', 'all'); %silencing warnings inside p618
    for i = 1:n_grid
        link_cfg = p618Config('Frequency', f, 'ElevationAngle', el_grid(i), ... 
            'Latitude', lat, 'Longitude', lon, 'TotalAnnualExceedance', 1, ... 
            'PolarizationTiltAngle', 45, 'AntennaDiameter', 0.5, 'AntennaEfficiency', 0.5);      
        [pl, ~, tsky] = p618PropagationLosses(link_cfg, 'StationHeight', 0);
        grid_Ag(i) = pl.Ag; grid_Ac(i) = pl.Ac; grid_Ar(i) = pl.Ar;
        grid_As(i) = pl.As; grid_At(i) = pl.At; grid_Tsky(i) = tsky;
    end
    warning(old_warn_state); % getting back at the old warning state
    
    loss_struct.Ag = interp1(el_grid, grid_Ag, el_vec, 'linear');
    loss_struct.Ac = interp1(el_grid, grid_Ac, el_vec, 'linear');
    loss_struct.Ar = interp1(el_grid, grid_Ar, el_vec, 'linear');
    loss_struct.As = interp1(el_grid, grid_As, el_vec, 'linear');
    loss_struct.At = interp1(el_grid, grid_At, el_vec, 'linear');
    sky_temp_K = interp1(el_grid, grid_Tsky, el_vec, 'linear');
end

function [loss_rx, loss_tx] = Steering_loss_calc(el_vec, range_vec, tx_type, rx_type)
    Re = 6378.14e3;     
    cos_exponent = 1.5; 
    theta_rx = 90 - el_vec; 
    R_sat_sq = Re^2 + range_vec.^2 + 2 * Re .* range_vec .* sind(el_vec);
    R_sat = sqrt(R_sat_sq);
    sin_theta_tx = (Re .* cosd(el_vec)) ./ R_sat;
    theta_tx = asind(sin_theta_tx); 
    
    loss_tx = 0; loss_rx = 0;
    if contains(tx_type, 'array')
        loss_tx = -10 * log10(cosd(theta_tx).^cos_exponent);
    end
    if contains(rx_type, 'array')
        loss_rx = -10 * log10(cosd(theta_rx).^cos_exponent);
    end
end

function N0_dBmHz = Noise_density_calc(antenna_temp, NF)
    k_dB = 10*log10(1.38e-23); 
    F_lin = 10^(NF/10);
    T_rx = (F_lin - 1) * 290;
    T_sys = antenna_temp + T_rx;
    N0_dBmHz = (10*log10(T_sys) + k_dB) + 30;
end


