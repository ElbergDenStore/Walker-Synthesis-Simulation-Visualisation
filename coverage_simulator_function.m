function metrics = coverage_simulator_function(Cfg, use_parallel, calc_link)
% RUN_SATELLITE_SIM Simulates satellite coverage and link budget.
    fprintf('\n Starting Simulation: %d Sats, %.1f deg Inclination\n', Cfg.Total_sats, Cfg.Inclination);
    tic

    %% Scenario & Constellation Setup
    sc = satelliteScenario;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;

    r_earth = 6378.14e3;
    % Coverage reference latitude for seam-ratio geometry
    if isfield(Cfg, 'Lat_range_deg')
        min_lat_cov = min(Cfg.Lat_range_deg);
    else
        min_lat_cov = 0;
    end
    fprintf("Constructing Satellites\n")
    if Cfg.WalkerStar == true
        sats = asymmetrical_walker_star_generation(sc, Cfg.Orbit_height, Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Min_elevation_UE, "two-body-keplerian", min_lat_cov);
    else
        sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="two-body-keplerian");
    end

    % Get the massive ECEF matrix instantly from SGP4
    fprintf("Propagating Satellites\n")
    [sat_pos_raw, ~, simTimes] = states(sats, "CoordinateFrame", "ECEF");
    
    % Permute to [3 x NumSats x nT] to make the implicit expansion math easy
    sat_pos_ecef = permute(sat_pos_raw, [1, 3, 2]);
    
    % Cache the dimensions for the math loop
    num_sats = size(sat_pos_ecef, 2);
    nT = size(sat_pos_ecef, 3);
    
    %% Create the UEs Struct Array (Pre-allocated)
    if ~isfield(Cfg, 'Flat_UE_array') || ~isfield(Cfg.Flat_UE_array, 'Lats') || ~isfield(Cfg.Flat_UE_array, 'Lons')
        error('Cfg.Flat_UE_array with fields Lats and Lons is required. Generate UEs before calling coverage_simulator_function.');
    end

    UE_lats = Cfg.Flat_UE_array.Lats;
    UE_lons = Cfg.Flat_UE_array.Lons;
    UE_lats = UE_lats(:);
    UE_lons = UE_lons(:);

    ue_pos_ecef = lla2ecef([UE_lats, UE_lons, zeros(length(UE_lats), 1)]);
    Cfg.NumUEs = length(UE_lats);
    
    % 1. Define the perfectly sized SimData template
    empty_SimData = struct(...
        'Time',          NaN(1, nT), ... 
        'SatID',         NaN(1, nT), ...
        'Range',         NaN(1, nT), ...
        'Elevation_deg', NaN(1, nT), ...
        'Azimuth_deg',   NaN(1, nT), ...
        'Num_visible',   zeros(1, nT) ...
    );

    % 2. Conditionally define the Link template
    if calc_link
        empty_Absorption = struct('Ag', NaN(1, nT), 'Ac', NaN(1, nT), ...
                                  'Ar', NaN(1, nT), 'As', NaN(1, nT), 'At', NaN(1, nT));
        empty_link = struct(...
            'Frequency',         NaN, ...
            'Bandwidth',         NaN, ...
            'FSPL',              NaN(1, nT), ...
            'Absorption',        empty_Absorption, ...
            'T_antenna',         NaN(1, nT), ...
            'Rx_steering_loss',  NaN(1, nT), ...
            'Tx_steering_loss',  NaN(1, nT), ...
            'Total_loss',        NaN(1, nT), ...
            'Adjusted_EIRP_density_dBmHz', NaN(1, nT), ...
            'PFD_W_MHz',         NaN(1, nT), ...
            'Noise_density_dBmHz', NaN(1, nT), ...
            'Carrier_density_dBmHz', NaN(1, nT), ...
            'Interference_density_dBmHz', NaN(1, nT), ...
            'P_noise',           NaN(1, nT), ...
            'Rx_Power',          NaN(1, nT), ...
            'SNR',               NaN(1, nT), ...
            'SIR',               NaN(1, nT), ...
            'SINR',              NaN(1, nT), ...
            'Throughput',        NaN(1, nT) ...
        );
    end
    
    % 3. Lock in memory for the Struct Array
    UEs(Cfg.NumUEs).Lat = []; 
    for idx = 1:Cfg.NumUEs
        UEs(idx).Lat = UE_lats(idx);
        UEs(idx).Lon = UE_lons(idx);
        UEs(idx).Name = sprintf('UE%d', idx);
        UEs(idx).SimData = empty_SimData;
        
        if calc_link
            if isfield(Cfg, 'DL')
                UEs(idx).DL = empty_link;
                UEs(idx).DL.Frequency = Cfg.DL.f;
                UEs(idx).DL.Bandwidth = Cfg.DL.B;
            end
            if isfield(Cfg, 'UL')
                UEs(idx).UL = empty_link;
                UEs(idx).UL.Frequency = Cfg.UL.f;
                UEs(idx).UL.Bandwidth = Cfg.UL.B;
            end
        end
    end
    
    %% Coverage Simulation
    dq = parallel.pool.DataQueue;
    updateLiveScriptProgress(Cfg.NumUEs, true); 
    afterEach(dq, @(~) updateLiveScriptProgress(Cfg.NumUEs, false));
    if use_parallel
        num_workers = Inf; 
    else
        num_workers = 0;   
    end
    
    min_elevation_UE = Cfg.Min_elevation_UE;

    Cfg.Num_workers = num_workers; % Will be used in link calculation
    
    tic
    parfor (idx = 1:Cfg.NumUEs, num_workers)
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

        % ue = groundStation(sc, UEs(idx).Lat, UEs(idx).Lon);
        % 
        % [az_mat, el_mat, r_mat, simTimes] = aer(ue, sats);
        % 

        valid_mask = el_mat >= min_elevation_UE;
        Num_visible = sum(valid_mask,1); 
        has_service = Num_visible > 0;
        
        r_temp = r_mat;
        r_temp(~valid_mask) = Inf; 
        [best_ranges, best_sat_idx] = min(r_temp, [], 1); 
        
        % Write directly to pre-allocated slice
        UEs(idx).SimData.Num_visible = Num_visible;
        UEs(idx).SimData.Time = simTimes;
        
        
        if any(has_service)
            UEs(idx).SimData.Range(has_service) = best_ranges(has_service);
            
            best_sats_valid = best_sat_idx(has_service);
            UEs(idx).SimData.SatID(has_service) = best_sats_valid;
            
            valid_cols = find(has_service);
            num_rows = size(el_mat, 1);
            lin_idxs = best_sats_valid + (valid_cols - 1) * num_rows;
            
            UEs(idx).SimData.Elevation_deg(has_service) = el_mat(lin_idxs);
            UEs(idx).SimData.Azimuth_deg(has_service)   = az_mat(lin_idxs);
        end
        
        send(dq, []);
    end
    fprintf('\nGeometry calculation complete (%.1f sec).\n', toc);
    
    %% Coverage Stats (Fully Vectorized)
    SimDataArray = [UEs.SimData];
    counts = vertcat(SimDataArray.Num_visible);
    
    more_1_satellites = (counts >= 1);
    prob_coverage = 100 * sum(more_1_satellites, 2) ./ nT; 
    minNumberSatellites = min(counts, [], 2);
    meanNumberSatellites = mean(counts, 2);
    
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
    metrics.Num_visible            = counts; 
    
    % Store the final struct array in metrics
    metrics.SimData = SimDataArray;
    
    % Fallbacks if Link Calc is off
    metrics.throughput_10pct  = NaN;
    metrics.throughput_mean   = NaN;
    
    %% Link Budget Calculation (2D Vectorized & Cell-Batched)
    if calc_link
        tic
        % 1. Extract massive 2D matrices from the struct array
        el_mat    = single(vertcat(SimDataArray.Elevation_deg)); 
        az_mat    = single(vertcat(SimDataArray.Azimuth_deg));   
        range_mat = single(vertcat(SimDataArray.Range));      
        lat_vec   = [UEs.Lat]';
        lon_vec   = [UEs.Lon]';
        
        if isfield(Cfg, 'DL')
            batch_size = 100;
            num_batches = ceil(Cfg.NumUEs / batch_size);
            
            % Pre-allocate a CELL ARRAY to safely catch parallel worker outputs
            batch_results = cell(num_batches, 1);
            
            fprintf('\nProcessing Link Budget in %d batches...\n', num_batches);
            
            % --- PARALLEL BATCH PROCESSING ---
            parfor (b = 1:num_batches, num_workers)
                start_idx = (b - 1) * batch_size + 1;
                end_idx   = min(b * batch_size, Cfg.NumUEs);
                
                % Slice the inputs for the current chunk
                el_chunk    = el_mat(start_idx:end_idx, :);
                az_chunk    = az_mat(start_idx:end_idx, :);
                range_chunk = range_mat(start_idx:end_idx, :);
                lat_chunk   = lat_vec(start_idx:end_idx);
                lon_chunk   = lon_vec(start_idx:end_idx);
                
                % Run the calculation and store it in the sliced cell array
                batch_results{b} = link_calc_matrix(el_chunk, az_chunk, range_chunk, lat_chunk, lon_chunk, Cfg.DL, Cfg);
            end
            
            % --- SEQUENTIAL UNPACKING ---
            % distribute the cell results back into the UEs struct
            for b = 1:num_batches
                start_idx = (b - 1) * batch_size + 1;
                end_idx   = min(b * batch_size, Cfg.NumUEs);
                chunk     = batch_results{b};
                
                local_idx = 1;
                for idx = start_idx:end_idx
                    UEs(idx).DL.Frequency         = chunk.Frequency;
                    UEs(idx).DL.Bandwidth         = chunk.Bandwidth;
                    UEs(idx).DL.FSPL              = chunk.FSPL(local_idx, :);
                    UEs(idx).DL.Total_loss        = chunk.Total_loss(local_idx, :);
                    UEs(idx).DL.PFD_W_MHz         = chunk.PFD_W_MHz(local_idx, :);
                    UEs(idx).DL.Noise_density_dBmHz = chunk.Noise_density_dBmHz(local_idx, :);
                    UEs(idx).DL.Carrier_density_dBmHz = chunk.Carrier_density_dBmHz(local_idx, :);
                    UEs(idx).DL.Interference_density_dBmHz = chunk.Interference_density_dBmHz(local_idx, :);
                    UEs(idx).DL.SNR               = chunk.SNR(local_idx, :);
                    UEs(idx).DL.SIR               = chunk.SIR(local_idx, :);
                    UEs(idx).DL.SINR              = chunk.SINR(local_idx, :);
                    UEs(idx).DL.serving_beam_idx     = chunk.serving_beam_idx(local_idx, :);
                    UEs(idx).DL.serving_beam_signal_lin = chunk.serving_beam_signal_lin(local_idx, :);
                    UEs(idx).DL.interference_lin   = chunk.interference_lin(local_idx, :);
                    UEs(idx).DL.Absorption.At     = chunk.Absorption_At(local_idx, :);
                    
                    local_idx = local_idx + 1;
                end
            end

            DL_structs = [UEs.DL];
            DL_SINR = vertcat(DL_structs.SINR);
            DL_serving_beam_idx = vertcat(DL_structs.serving_beam_idx);
            SimData_structs = [UEs.SimData];
            SimData_SatID = vertcat(SimData_structs.SatID);
            DL_Throughput = calculate_throughput_matrix(DL_SINR, DL_serving_beam_idx, SimData_SatID, Cfg.DL.BeamGrid, Cfg.DL.B, Cfg.Share_bandwidth, Cfg.Modified_shannon);

            for idx = 1:Cfg.NumUEs %% TODO look closely at this, it looks stupid and inefficient, messy
                UEs(idx).DL.Throughput = DL_Throughput(idx, :);
            end
        end
        fprintf('\nLoss calculation complete (%.1f sec).\n', toc);
    end
    
    %% ====================================================================
    %% PACKAGE FINAL METRICS (The Pure Function Output)
    %% ====================================================================
    
    metrics.worst_coverage_percent = min(prob_coverage);
    metrics.Num_visible            = counts; 
    metrics.minNumberSatellites    = minNumberSatellites;
    metrics.meanNumberSatellites   = meanNumberSatellites;
    metrics.prob_coverage          = prob_coverage;
    metrics.UEs                    = UEs; % Contains SimData and DL/UL links
    metrics.Cfg                    = Cfg;
    
end