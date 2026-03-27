function metrics = fast_coverage_simulator_function(Cfg, reset_cache, use_parallel, calc_link)
% FAST_COVERAGE_SIMULATOR_FUNCTION Optmized simulator strictly using ECEF math.
    fprintf('\n Starting Fast Simulation: %d Sats, %.1f deg Inclination\n', Cfg.Total_sats, Cfg.Inclination);
    tic
    %% 2. Constellation Setup (Cached)
    persistent cached_sc cached_UE_lats cached_UE_lons cached_ue_pos_ecef last_Cfg
    
    if reset_cache
        cached_sc = [];
        cached_UE_lats = [];
        cached_UE_lons = [];
        cached_ue_pos_ecef = [];
        last_Cfg = [];
        fprintf(' [!] Cache manually reset for fresh initialization.\n');
    end
    
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
    
    
    % Reshape Satellites to [3 x NumSats x nT x 1]
    sat_pos_4d = reshape(sat_pos_ecef, 3, num_sats, nT, 1);
    
    % Reshape UEs to [3 x 1 x 1 x NumUEs]
    ue_pos_4d = reshape(ue_pos_ecef', 3, 1, 1, Cfg.NumUEs);
    
    % Subtract instantly! Result is [3 x NumSats x nT x NumUEs]
    vec_ecef_4d = sat_pos_4d - ue_pos_4d;
    
    % ===============================================================
    % 2. 4D ROTATION (ECEF TO ENU) via pagemtimes
    % ===============================================================
    % Calculate sines and cosines for all UEs simultaneously [1 x NumUEs]
    slat = sind(UE_lats'); clat = cosd(UE_lats');
    slon = sind(UE_lons'); clon = cosd(UE_lons');
    
    % Build a 3D array of rotation matrices [3 x 3 x NumUEs]
    R_ecef_to_enu = zeros(3, 3, Cfg.NumUEs);
    R_ecef_to_enu(1,1,:) = -slon;           R_ecef_to_enu(1,2,:) = clon;            R_ecef_to_enu(1,3,:) = 0;
    R_ecef_to_enu(2,1,:) = -slat.*clon;     R_ecef_to_enu(2,2,:) = -slat.*slon;     R_ecef_to_enu(2,3,:) = clat;
    R_ecef_to_enu(3,1,:) = clat.*clon;      R_ecef_to_enu(3,2,:) = clat.*slon;      R_ecef_to_enu(3,3,:) = slat;
    
    % Reshape vec_ecef so we can multiply it: [3 x (NumSats*nT) x NumUEs]
    vec_ecef_pages = reshape(vec_ecef_4d, 3, num_sats * nT, Cfg.NumUEs);
    
    % Multiply all 12 rotation matrices by their respective coordinate blocks instantly!
    vec_enu_pages = pagemtimes(R_ecef_to_enu, vec_ecef_pages);
    
    % Reshape back to our beautiful 4D Tensor: [3 x NumSats x nT x NumUEs]
    vec_enu_4d = reshape(vec_enu_pages, 3, num_sats, nT, Cfg.NumUEs);
    
    % ===============================================================
    % 3. EXTRACT AER (Vectorized across all 4 Dimensions)
    % ===============================================================
    E = vec_enu_4d(1, :, :, :);
    N = vec_enu_4d(2, :, :, :);
    U = vec_enu_4d(3, :, :, :);
    
    r_4d = sqrt(E.^2 + N.^2 + U.^2);
    el_4d = asind(U ./ r_4d);
    
    az_4d = atan2d(E, N);
    az_4d(az_4d < 0) = az_4d(az_4d < 0) + 360;
    
    % Squeeze down to 3D matrices: [NumSats x nT x NumUEs]
    r_mat_all = reshape(r_4d, num_sats, nT, Cfg.NumUEs);
    el_mat_all = reshape(el_4d, num_sats, nT, Cfg.NumUEs);
    az_mat_all = reshape(az_4d, num_sats, nT, Cfg.NumUEs);

    % ===============================================================
    % 4. ASSIGN TO STRUCT (Fast slicing loop)
    % ===============================================================
    % We still use a loop to pack the struct array, but NO math happens here.
    min_el = Cfg.Min_elevation_UE;
        
   for idx = 1:Cfg.NumUEs
        % Slice out this UE's completed matrices
        r_mat  = r_mat_all(:, :, idx);
        el_mat = el_mat_all(:, :, idx);
        az_mat = az_mat_all(:, :, idx);
        
        valid_mask = el_mat >= min_el;
        Num_visible = sum(valid_mask, 1); 
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
