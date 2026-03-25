function Cfg = get_cfg(height_km, constellation_type, ue_grid_size, duration, frequency)
    arguments %Set default behaviors if not all inputs are input
        height_km (1,1) double = 1000
        constellation_type (1,1) string = "walkerdelta"
        ue_grid_size (1,1) string = "small"
        duration (1,1) string = "short"
        frequency (1,1) string = "ku"
    end

    %%%%% CONSTELLATION %%%%
    optimal_constellation = load("optimal_constellations.mat");

    constellation_idx = find(optimal_constellation.heights_km >= height_km, 1, 'first');
    switch lower(constellation_type)
        case 'walkerdelta'
            Cfg = optimal_constellation.best_delta_sats(constellation_idx);
            Cfg.WalkerStar = false;
        case 'walkerstar'
            Cfg = optimal_constellation.star_sats(constellation_idx);
            Cfg.WalkerStar = true;
        otherwise
            error('Invalid constellation type');
    end
    Cfg.Total_sats = Cfg.Sats_per_plane * Cfg.Num_planes;
    Cfg.SampleTime = 60; % seconds
    Cfg.Min_elevation_UE = 20;
    Cfg.Orbit_height = height_km*1e3;


    %%%%% UE GRID SIZE %%%%
    Cfg.Lon_vec       = linspace(-60, 30, 2);
    Cfg.Equal_UE_area = true;
    switch lower(ue_grid_size)
        case 'small'
            Cfg.Lat_vec = linspace(55, 85, 2); % 6
        case 'medium'
            Cfg.Lat_vec = linspace(55, 85, 4); % 20
        case 'big'
            Cfg.Lat_vec = linspace(55, 85, 8); % 73
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
            Cfg.DL.B         = 2e6;     
            Cfg.DL.f         = 20e9;
            Cfg.DL.Tx_type   = "array";      
            Cfg.DL.G_rx      = 36;      
            Cfg.DL.Rx_type   = "array";      
            Cfg.DL.NF        = 5;
        case 'ku'
            Cfg.Target_PFD_MHz = -115;
            Cfg.DL.B         = 2e6;     
            Cfg.DL.f         = 12e9;
            Cfg.DL.Tx_type   = "array";      
            Cfg.DL.G_rx      = 33;      
            Cfg.DL.Rx_type   = "array";      
            Cfg.DL.NF        = 5;
        otherwise
            error('Invalid frequency');
    end
    Cfg.DL.Direction = "DL";
    Cfg.DL.G_tx          = get_adjusted_tx_gain(Cfg.Orbit_height, Cfg.Min_elevation_UE, Cfg.DL.f); 
    Cfg.DL.Max_P_tx_dBm  = PFD_calc(Cfg.Target_PFD_MHz, Cfg.DL.G_tx, Cfg.DL.B, Cfg.Orbit_height, Cfg.Min_elevation_UE);
    Cfg.DL.Max_EIRP_dBm  = Cfg.DL.Max_P_tx_dBm + Cfg.DL.G_tx;

end