function metrics = fast_coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link)
% FAST_COVERAGE_SIMULATOR_FUNCTION Optmized simulator strictly using ECEF math.
    fprintf('\n Starting Fast Simulation: %d Sats, %.1f deg Inclination\n', Cfg.Total_sats, Cfg.Inclination);
    
    %% 2. Constellation Setup (Cached)
    persistent cached_sc
    
    if isempty(cached_sc) || ~isvalid(cached_sc)
        cached_sc = satelliteScenario;
    else
        if ~isempty(cached_sc.Satellites)
            delete(cached_sc.Satellites);
        end
        if ~isempty(cached_sc.GroundStations)
            delete(cached_sc.GroundStations);
        end
    end
    
    sc = cached_sc;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;
    r_earth = 6378.14e3;
    
    if Cfg.WalkerStar == true
        sats = asymmetrical_walker_star_generation(Cfg.Orbit_height, Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane);
    else
        sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");
    end
    
    % Get the massive ECEF matrix instantly from SGP4
    [sat_pos_raw, ~, simTimes] = states(sats, "CoordinateFrame", "ECEF");
    
    % Permute to [3 x NumSats x nT] to make our implicit expansion math easy
    sat_pos_ecef = permute(sat_pos_raw, [1, 3, 2]);
    
    % Cache the dimensions for the math loop
    num_sats = size(sat_pos_ecef, 2);
    nT = size(sat_pos_ecef, 3);
    
    %% 3. Create UEs Struct Array
    persistent cached_UE_lats cached_UE_lons cached_ue_pos_ecef last_Cfg
    
    recon_UEs = false;
    if isempty(cached_UE_lats)
        recon_UEs = true;
    elseif ~isequal(last_Cfg.Lat_vec, Cfg.Lat_vec) || ~isequal(last_Cfg.Lon_vec, Cfg.Lon_vec)
        recon_UEs = true;
    elseif last_Cfg.Equal_UE_area ~= Cfg.Equal_UE_area
        recon_UEs = true;
    end
    
    if recon_UEs
        if Cfg.Equal_UE_area == true
            [UE_lats, UE_lons] = generate_equal_ish_area_UEs(Cfg.Lat_vec, Cfg.Lon_vec);
        elseif Cfg.Accept_Flat_UE_array == true
            UE_lats = Cfg.Flat_UE_array.Lats;
            UE_lons = Cfg.Flat_UE_array.Lons;
        else
            [UE_lats, UE_lons] = meshgrid(Cfg.Lat_vec, Cfg.Lon_vec);
        end
        cached_UE_lats = UE_lats(:);
        cached_UE_lons = UE_lons(:);
        cached_ue_pos_ecef = lla2ecef([cached_UE_lats, cached_UE_lons, zeros(length(cached_UE_lats), 1)]);
        
        last_Cfg.Lat_vec = Cfg.Lat_vec;
        last_Cfg.Lon_vec = Cfg.Lon_vec;
        last_Cfg.Equal_UE_area = Cfg.Equal_UE_area;
    end
    
    UE_lats = cached_UE_lats;
    UE_lons = cached_UE_lons;
    ue_pos_ecef = cached_ue_pos_ecef;
    Cfg.NumUEs = length(UE_lats);
    
    % Lock in memory for Struct Array
    empty_SimData = struct('Time', simTimes, 'SatID', NaN(1, nT), ...
        'Range', NaN(1, nT), 'Elevation_deg', NaN(1, nT), ...
        'Azimuth_deg', NaN(1, nT), 'Num_visible', zeros(1, nT));
        
    if calc_link
        empty_Absorption = struct('Ag', NaN(1, nT), 'Ac', NaN(1, nT), ...
                                  'Ar', NaN(1, nT), 'As', NaN(1, nT), 'At', NaN(1, nT));
        empty_link = struct('Frequency', NaN, 'Bandwidth', NaN, 'FSPL', NaN(1, nT), ...
            'Absorption', empty_Absorption, 'T_antenna', NaN(1, nT), ...
            'Rx_steering_loss', NaN(1, nT), 'Tx_steering_loss', NaN(1, nT), ...
            'Total_loss', NaN(1, nT), 'Adjusted_EIRP_dBm', NaN(1, nT), ...
            'PFD_W_MHz', NaN(1, nT), 'P_noise', NaN(1, nT), 'Rx_Power', NaN(1, nT), ...
            'SNR', NaN(1, nT), 'SIR', NaN(1, nT), 'SINR', NaN(1, nT), 'Throughput', NaN(1, nT));
    end
    
    UEs(Cfg.NumUEs).Lat = []; 
    for idx = 1:Cfg.NumUEs
        UEs(idx).Lat  = UE_lats(idx);
        UEs(idx).Lon  = UE_lons(idx);
        UEs(idx).Name = sprintf('UE%d', idx);
        UEs(idx).SimData = empty_SimData;
        
        if calc_link
            if isfield(Cfg, 'DL')
                UEs(idx).DL = empty_link;
                UEs(idx).DL.Frequency = Cfg.DL.f; UEs(idx).DL.Bandwidth = Cfg.DL.B;
            end
            if isfield(Cfg, 'UL')
                UEs(idx).UL = empty_link;
                UEs(idx).UL.Frequency = Cfg.UL.f; UEs(idx).UL.Bandwidth = Cfg.UL.B;
            end
        end
    end
    
    %% Fast Coverage Simulation using pure ECEF math
    % if use_parallel, num_workers = Inf; else, num_workers = 0; end
    min_elevation_UE = Cfg.Min_elevation_UE;
    % Cfg.Num_workers = num_workers;
    
    tic
    for idx = 1:Cfg.NumUEs
        ue_xyz = ue_pos_ecef(idx, :)'; % 3x1 vector
        lat = UE_lats(idx);
        lon = UE_lons(idx);
        
        % 2. Vector from UE to ALL satellites at ALL times [3 x NumSats x nT]
        vec_ecef = sat_pos_ecef - ue_xyz;
        
        % ===============================================================
        % 3. ECEF TO ENU ROTATION (To get perfect Azimuth and Elevation)
        % ===============================================================
        slat = sind(lat); clat = cosd(lat);
        slon = sind(lon); clon = cosd(lon);
        
        % Standard transformation matrix from Earth-Centered to Local East-North-Up
        R_ecef_to_enu = [
            -slon,           clon,          0;
            -slat*clon,     -slat*slon,     clat;
             clat*clon,      clat*slon,     slat
        ];
        
        % Flatten the vectors, rotate them all instantly, and re-fold the matrix
        vec_ecef_flat = reshape(vec_ecef, 3, []);
        vec_enu_flat = R_ecef_to_enu * vec_ecef_flat;
        vec_enu = reshape(vec_enu_flat, 3, num_sats, nT);
        
        % Extract East, North, and Up components (Safely reshaping to preserve dimensions)
        E = reshape(vec_enu(1, :, :), num_sats, nT);
        N = reshape(vec_enu(2, :, :), num_sats, nT);
        U = reshape(vec_enu(3, :, :), num_sats, nT);
        
        % ===============================================================
        % 4. EXTRACT AER (Azimuth, Elevation, Range)
        % ===============================================================
        r_mat = sqrt(E.^2 + N.^2 + U.^2);
        el_mat = asind(U ./ r_mat);
        
        az_mat = atan2d(E, N);
        az_mat(az_mat < 0) = az_mat(az_mat < 0) + 360;
        
        valid_mask = el_mat >= min_elevation_UE;
        Num_visible = sum(valid_mask,1); 
        has_service = Num_visible > 0;
        
        r_temp = r_mat;
        r_temp(~valid_mask) = Inf; 
        [best_ranges, best_sat_idx] = min(r_temp, [], 1); 
        
        UEs(idx).SimData.Num_visible = Num_visible;
        
        if any(has_service)
            UEs(idx).SimData.Range(has_service) = best_ranges(has_service);
            
            best_sats_valid = best_sat_idx(has_service);
            UEs(idx).SimData.SatID(has_service) = best_sats_valid;
            
            valid_cols = find(has_service);
            lin_idxs = best_sats_valid + (valid_cols - 1) * num_sats;
            
            UEs(idx).SimData.Elevation_deg(has_service) = el_mat(lin_idxs);
            UEs(idx).SimData.Azimuth_deg(has_service)   = az_mat(lin_idxs);
        end
    end
    fprintf('\nGeometry calculation complete (%.1f sec).\n', toc);
    
    %% Coverage Stats (Fully Vectorized)
    SimDataArray = [UEs.SimData];
    counts = vertcat(SimDataArray.Num_visible);
    
    more_1_satellites = (counts >= 1);
    prob_coverage = 100 * sum(more_1_satellites, 2) ./ nT; 
    minNumberSatellites = min(counts, [], 2);
    meanNumberSatellites = mean(counts, 2);
    
    maxGapMinutes = zeros(Cfg.NumUEs,1);
    for idx = 1:Cfg.NumUEs
        row = (counts(idx,:) == 0);  
        maxZeroStreak = 0; currentStreak = 0;
        for k = 1:length(row)
            if row(k)
                currentStreak = currentStreak + 1;
                maxZeroStreak = max(maxZeroStreak, currentStreak);
            else
                currentStreak = 0;
            end
        end
        maxGapMinutes(idx) = maxZeroStreak * Cfg.SampleTime / 60;
    end
    
    metrics.worst_coverage_percent = min(prob_coverage);
    metrics.worst_gap_minutes      = max(maxGapMinutes);
    metrics.Num_visible            = counts; 
    metrics.SimData                = SimDataArray;
    
    metrics.throughput_10pct  = NaN;
    metrics.throughput_mean   = NaN;
    meanThroughput = NaN(Cfg.NumUEs, 1);
    
    %% Link Budget Calculation
    if calc_link
        tic
        el_mat_calc = vertcat(SimDataArray.Elevation_deg); 
        az_mat_calc = vertcat(SimDataArray.Azimuth_deg);   
        range_mat_calc = vertcat(SimDataArray.Range);      
        lat_vec_calc = [UEs.Lat]';
        lon_vec_calc = [UEs.Lon]';
        
        if isfield(Cfg, 'DL')
            DL_Result = link_calc_matrix(el_mat_calc, az_mat_calc, range_mat_calc, lat_vec_calc, lon_vec_calc, Cfg.DL, Cfg);
            for idx = 1:Cfg.NumUEs
                UEs(idx).DL.Frequency         = DL_Result.Frequency;
                UEs(idx).DL.Bandwidth         = DL_Result.Bandwidth;
                UEs(idx).DL.FSPL              = DL_Result.FSPL(idx, :);
                UEs(idx).DL.Total_loss        = DL_Result.Total_loss(idx, :);
                UEs(idx).DL.Adjusted_EIRP_dBm = DL_Result.Adjusted_EIRP_dBm(idx, :);
                UEs(idx).DL.PFD_W_MHz         = DL_Result.PFD_W_MHz(idx, :);
                UEs(idx).DL.SNR               = DL_Result.SNR(idx, :);
                UEs(idx).DL.Throughput        = DL_Result.Throughput(idx, :);
                UEs(idx).DL.Absorption.At     = DL_Result.Absorption_At(idx, :);
            end
            
            valid_thpt = DL_Result.Throughput(~isnan(DL_Result.Throughput));
            metrics.throughput_10pct = prctile(valid_thpt, 10);
            metrics.throughput_mean  = mean(valid_thpt);
        end
        fprintf('\nLoss calculation complete (%.1f sec).\n', toc);
    end
end
