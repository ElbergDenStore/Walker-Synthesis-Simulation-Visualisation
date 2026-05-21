function metrics = fast_coverage_simulator_function(Cfg, reset_cache, calc_link, use_SGP)
% FAST_COVERAGE_SIMULATOR_FUNCTION Optmized simulator strictly using ECEF math.
    fprintf('\n Starting Fast Simulation: %d Sats, %.1f deg Inclination\n', Cfg.Total_sats, Cfg.Inclination);
    tic

    % Coverage reference latitude for seam-ratio geometry
    if isfield(Cfg, 'Lat_range_deg')
        min_lat_cov = min(Cfg.Lat_range_deg);
    else
        min_lat_cov = 0;
    end

    %% 2. Constellation & Time Setup
    persistent cached_sc last_Cfg
    
    if reset_cache
        cached_sc = [];
        last_Cfg = [];
        fprintf(' [!] Cache manually reset for fresh initialization.\n');
    end

    if use_SGP
        %% SGP4 PROPAGATOR (Standard MATLAB Toolbox)
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
            sats = asymmetrical_walker_star_generation(sc, Cfg.Orbit_height, Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Min_elevation_UE, "sgp4", min_lat_cov);
        else
            sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
                Cfg.Inclination, Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
                Name="S4D", OrbitPropagator="sgp4");
        end
        
        % Get the massive ECEF matrix instantly from SGP4
        [sat_pos_raw, ~, simTimes] = states(sats, "CoordinateFrame", "ECEF");
        
        % Permute to [3 x NumSats x nT] to make our implicit expansion math easy
        sat_pos_ecef = permute(sat_pos_raw, [1, 3, 2]);
    else
        %% PURE MATH WALKER GENERATOR
        if Cfg.WalkerStar == true
            % Pure-math generator only supports Walker Delta; fall back to SGP4.
            fprintf(' [!] WalkerStar requested in pure-math mode – falling back to SGP4.\n');
            if isempty(cached_sc) || ~isvalid(cached_sc)
                cached_sc = satelliteScenario;
            else
                if ~isempty(cached_sc.Satellites),     delete(cached_sc.Satellites);     end
                if ~isempty(cached_sc.GroundStations), delete(cached_sc.GroundStations); end
            end
            sc = cached_sc;
            sc.StartTime  = Cfg.StartTime;
            sc.StopTime   = Cfg.StopTime;
            sc.SampleTime = Cfg.SampleTime;
            sats = asymmetrical_walker_star_generation(sc, Cfg.Orbit_height, Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Min_elevation_UE, "sgp4", min_lat_cov);
            [sat_pos_raw, ~, simTimes] = states(sats, "CoordinateFrame", "ECEF");
            sat_pos_ecef = permute(sat_pos_raw, [1, 3, 2]);
        else
            % Build the exact time vector using double math
            total_duration_sec = seconds(Cfg.StopTime - Cfg.StartTime);
            time_steps_sec = 0 : Cfg.SampleTime : total_duration_sec;

            simTimes = Cfg.StartTime + seconds(time_steps_sec);
            simTimes.TimeZone = 'UTC';

            sat_pos_ecef = fast_walker_ecef(Cfg.Orbit_height, Cfg.Inclination, ...
                                            Cfg.Num_planes, Cfg.Sats_per_plane, ...
                                            Cfg.Phasing, time_steps_sec, Cfg.StartTime);
        end
    end
    
    % Cache the dimensions for the math loop
    num_sats = size(sat_pos_ecef, 2);
    nT = size(sat_pos_ecef, 3);
    
    %% 3. Create UEs Struct Array

    if ~isfield(Cfg, 'Flat_UE_array') || ~isfield(Cfg.Flat_UE_array, 'Lats') || ~isfield(Cfg.Flat_UE_array, 'Lons')
        error('Cfg.Flat_UE_array with fields Lats and Lons is required. Generate UEs before calling coverage_simulator_function.');
    end
    UE_lats = Cfg.Flat_UE_array.Lats(:);
    UE_lons = Cfg.Flat_UE_array.Lons(:);
    ue_pos_ecef = lla2ecef([UE_lats, UE_lons, zeros(length(UE_lats), 1)]);
    Cfg.NumUEs = size(ue_pos_ecef, 1);
    
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
    % minNumberSatellites = min(counts, [], 2);
    % meanNumberSatellites = mean(counts, 2);
    
    % maxGapMinutes = zeros(Cfg.NumUEs,1);
    % for idx = 1:Cfg.NumUEs
    %     row = (counts(idx,:) == 0);  
    %     maxZeroStreak = 0; currentStreak = 0;
    %     for k = 1:length(row)
    %         if row(k)
    %             currentStreak = currentStreak + 1;
    %             maxZeroStreak = max(maxZeroStreak, currentStreak);
    %         else
    %             currentStreak = 0;
    %         end
    %     end
    %     maxGapMinutes(idx) = maxZeroStreak * Cfg.SampleTime / 60;
    % end
    
    metrics.worst_coverage_percent = min(prob_coverage);
    % metrics.worst_gap_minutes      = max(maxGapMinutes);
    % metrics.Num_visible            = counts; 
    % metrics.SimData                = SimDataArray;
    
    metrics.throughput_10pct  = NaN;
    metrics.throughput_mean   = NaN;
    meanThroughput = NaN(Cfg.NumUEs, 1);
    
    % %% Link Budget Calculation
    % if calc_link
    %     tic
    %     el_mat_calc = vertcat(SimDataArray.Elevation_deg); 
    %     az_mat_calc = vertcat(SimDataArray.Azimuth_deg);   
    %     range_mat_calc = vertcat(SimDataArray.Range);      
    %     lat_vec_calc = [UEs.Lat]';
    %     lon_vec_calc = [UEs.Lon]';
        
    %     if isfield(Cfg, 'DL')
    %         DL_Result = link_calc_matrix(el_mat_calc, az_mat_calc, range_mat_calc, lat_vec_calc, lon_vec_calc, Cfg.DL, Cfg);
    %         for idx = 1:Cfg.NumUEs
    %             UEs(idx).DL.Frequency         = DL_Result.Frequency;
    %             UEs(idx).DL.Bandwidth         = DL_Result.Bandwidth;
    %             UEs(idx).DL.FSPL              = DL_Result.FSPL(idx, :);
    %             UEs(idx).DL.Total_loss        = DL_Result.Total_loss(idx, :);
    %             UEs(idx).DL.Adjusted_EIRP_dBm = DL_Result.Adjusted_EIRP_dBm(idx, :);
    %             UEs(idx).DL.PFD_W_MHz         = DL_Result.PFD_W_MHz(idx, :);
    %             UEs(idx).DL.SNR               = DL_Result.SNR(idx, :);
    %             UEs(idx).DL.Throughput        = DL_Result.Throughput(idx, :);
    %             UEs(idx).DL.Absorption.At     = DL_Result.Absorption_At(idx, :);
    %         end
            
    %         valid_thpt = DL_Result.Throughput(~isnan(DL_Result.Throughput));
    %         metrics.throughput_10pct = prctile(valid_thpt, 10);
    %         metrics.throughput_mean  = mean(valid_thpt);
    %     end
    %     fprintf('\nLoss calculation complete (%.1f sec).\n', toc);
    % end
end

% =========================================================================
% PURE MATH WALKER GENERATOR
% =========================================================================
function sat_pos_ecef = fast_walker_ecef(Orbit_height, Inc_deg, P, S, F, time_steps_sec, StartTime)
    % Calculates the exact ECEF coordinates of a Walker Delta constellation
    
    T = P * S;
    a = 6378.137e3 + Orbit_height; % Semi-major axis (meters)
    mu = 3.986004418e14;           % Earth's gravitational constant
    we = 7.2921150e-5;             % Earth's rotation rate (rad/s)
    inc = deg2rad(Inc_deg);
    
    n = sqrt(mu / a^3);            % Mean motion (rad/s)
    nT = length(time_steps_sec);
    
    % --- THE MISSING LINK: EARTH'S INITIAL ROTATION (GMST) ---
    % 1. Convert StartTime to Julian Date
    JD = juliandate(StartTime);
    D = JD - 2451545.0; % Days since Jan 1, 2000, 12:00 UTC
    
    % 2. Calculate Greenwich Mean Sidereal Time in degrees
    GMST_deg = mod(280.46061837 + 360.98564736629 * D, 360);
    theta_g0 = deg2rad(GMST_deg);
    
    % Preallocate the [3 x NumSats x nT] matrix
    sat_pos_ecef = zeros(3, T, nT);
    
    sat_idx = 1;
    for p = 0:(P-1)
        RAAN = p * (2*pi / P); % Right Ascension of the Ascending Node
        
        for s = 0:(S-1)
            % Initial Phase (Mean Anomaly)
            M0 = s * (2*pi / S) + p * F * (2*pi / T);
            
            % Angle within the orbital plane over time (1 x nT vector)
            theta = M0 + n * time_steps_sec; 
            
            % 1. Position in the 2D orbital plane
            x_orb = a * cos(theta);
            y_orb = a * sin(theta);
            
            % 2. Rotate to 3D Earth-Centered Inertial (ECI)
            X_eci = x_orb * cos(RAAN) - y_orb * cos(inc) * sin(RAAN);
            Y_eci = x_orb * sin(RAAN) + y_orb * cos(inc) * cos(RAAN);
            Z_eci = y_orb * sin(inc);
            
            % 3. Rotate to Earth-Centered Earth-Fixed (ECEF) by spinning the Earth!
            % We add the initial offset (theta_g0) to the rotation over time
            theta_g = theta_g0 + we * time_steps_sec; 
            
            X_ecef = X_eci .* cos(theta_g) + Y_eci .* sin(theta_g);
            Y_ecef = -X_eci .* sin(theta_g) + Y_eci .* cos(theta_g);
            Z_ecef = Z_eci; % Z (North Pole) is unaffected by rotation
            
            % Store in our 3D tensor
            sat_pos_ecef(1, sat_idx, :) = X_ecef;
            sat_pos_ecef(2, sat_idx, :) = Y_ecef;
            sat_pos_ecef(3, sat_idx, :) = Z_ecef;
            
            sat_idx = sat_idx + 1;
        end
    end
end
