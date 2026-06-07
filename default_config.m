function Cfg = default_config(height_km, constellation_type, ue_grid_size, duration, frequency)
% DEFAULT_CONFIG  Build the default simulation configuration struct (Cfg).
%   Returns a fully-populated Cfg with constellation, link-budget, beam and
%   timing settings derived from HEIGHT_KM, CONSTELLATION_TYPE ("walkerdelta" or
%   "walkerstar"), UE_GRID_SIZE, DURATION and FREQUENCY band.  All arguments are
%   optional and default to a small Ku-band Walker Delta run.
    arguments %Set default behaviors if not all inputs are input
        height_km (1,1) double = 1000
        constellation_type (1,1) string = "walkerdelta"
        ue_grid_size (1,1) string = "small"
        duration (1,1) string = "short"
        frequency (1,1) string = "ku"
    end

    %%%%% CONSTELLATION %%%%
    Lat_range_deg = [54+(35/60), 83+(40/60)];
    Lon_range_deg = [-60, 30];
    min_elevation_UE = 20;
    requested_constellation = lower(constellation_type);
    mat_path = "optimal_constellations.mat";
    if isfile(mat_path)
        optimal_constellation = load(mat_path);
        [~, constellation_idx] = min(abs(optimal_constellation.heights_km - height_km));

        switch requested_constellation
            case 'walkerdelta'
                Cfg = optimal_constellation.best_delta_sats(constellation_idx, :);
                Cfg = table2struct(Cfg);
                Cfg.WalkerStar = false;
            case 'walkerstar'
                Cfg = optimal_constellation.star_sats(constellation_idx);
                Cfg.WalkerStar = true;
            otherwise
                error('Invalid constellation type');
        end
    else
        warning('No optimal_constellations.mat found. Switching to analytical Walker Star solution.');
        [Num_planes, Sats_per_plane, Total_sats] = calculate_walker_star(height_km, min(Lat_range_deg), min_elevation_UE);
        Cfg = struct('Orbit_height', height_km, ...
                     'Total_sats', Total_sats, ...
                     'Num_planes', Num_planes, ...
                     'Phasing', Num_planes / 2, ...
                     'Inclination', 87, ...
                     'Sats_per_plane', Sats_per_plane, ...
                     'WalkerStar', true);
    end

    Cfg.Total_sats = Cfg.Sats_per_plane * Cfg.Num_planes;
    Cfg.SampleTime = 60; % seconds
    Cfg.Min_elevation_UE = min_elevation_UE;
    Cfg.Lat_range_deg = Lat_range_deg;
    Cfg.Orbit_height = height_km*1e3;
    Cfg.FRF = 3;
    Cfg.RU = 1;
    Cfg.Use_P618 = false;
    Cfg.Modified_shannon = true; % Set to true for User/L2 achievable throughput (Modified Shannon), false for raw PHY Shannon capacity
    Cfg.Simple_Atmospheric_Loss_dB = 1;
    


    %%%%% UE GRID SIZE %%%%
    switch lower(ue_grid_size)
        case 'small'
            [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = generate_equal_area_ues(Lat_range_deg, Lon_range_deg, 6);
        case 'medium'
            [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = generate_equal_area_ues(Lat_range_deg, Lon_range_deg, 69);
        case 'big'
            [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = generate_equal_area_ues(Lat_range_deg, Lon_range_deg, 420);
        otherwise
            error('Invalid ue_grid_size');
    end

    %%%%% DURATION %%%%
    Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
    switch lower(duration)
        case 'short'
            Cfg.StopTime   = datetime('1-Jun-2025 12:59:59', 'TimeZone', 'UTC'); % 1 hour
        case 'medium'
            Cfg.StopTime   = datetime('1-Jun-2025 23:59:59', 'TimeZone', 'UTC'); % 12 hours
        case 'long'
            Cfg.StopTime   = datetime('3-Jun-2025 11:59:59', 'TimeZone', 'UTC'); % 48 hours

        otherwise
            error('Invalid duration');
    end

    %%%%% FREQUENCY %%%%
    Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
    switch lower(frequency)
        case 'fr1'
            % PFD_regulation = -113;
            Cfg.Target_PFD_MHz = -116;
            Cfg.DL.B         = 180e3;
            Cfg.DL.f         = 2.2e9;
            Cfg.DL.Tx_type   = "array";      
            Cfg.DL.G_rx      = 3; % Dipole antenna typically
            Cfg.DL.Rx_type   = "array";      
            Cfg.DL.NF        = 7; %  TR 38821
        case 'ka'
            Cfg.Target_PFD_MHz = -115;
            Cfg.DL.B         = 50e6;     
            Cfg.DL.f         = 20e9;
            Cfg.DL.Tx_type   = "array";      
            Cfg.DL.G_rx      = 36;      
            Cfg.DL.Rx_type   = "array";      
            Cfg.DL.NF        = 5;
        case 'ku'
            Cfg.Target_PFD_MHz = -115;
            Cfg.DL.B         = 50e6;     
            Cfg.DL.f         = 12e9;
            Cfg.DL.Tx_type   = "array";      
            Cfg.DL.G_rx      = 33;      
            Cfg.DL.Rx_type   = "array";      
            Cfg.DL.NF        = 5;
        otherwise
            error('Invalid frequency');
    end
    Cfg.DL.Direction     = "DL";
    Cfg.Share_bandwidth  = true; % All the bandwidth is shared per beam
    Cfg.DL.G_tx          = get_adjusted_tx_gain(Cfg.Orbit_height, Cfg.Min_elevation_UE, Cfg.DL.f); 
    Cfg.DL.Max_P_tx_dBm  = PFD_calc(Cfg.Target_PFD_MHz, Cfg.DL.G_tx, Cfg.DL.B, Cfg.Orbit_height, Cfg.Min_elevation_UE);
    Cfg.DL.Max_EIRP_dBm  = Cfg.DL.Max_P_tx_dBm + Cfg.DL.G_tx;
    Cfg.DL.Max_EIRP_dBm_Hz = Cfg.DL.Max_EIRP_dBm - 10*log10(Cfg.DL.B);

    % Configure generalized BeamGrid for standard constellations
    Cfg.DL.BeamGrid = calculate_hexagonal_beams(Cfg.DL.G_tx, Cfg.DL.f, Cfg.Orbit_height, Cfg.Min_elevation_UE, Cfg.DL.Max_EIRP_dBm_Hz, Cfg.FRF);

end