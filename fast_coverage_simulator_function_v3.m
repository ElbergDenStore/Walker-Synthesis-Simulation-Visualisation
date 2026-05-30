function metrics = fast_coverage_simulator_function_v3(Cfg, reset_cache, calc_link, use_SGP, cancel_queue)
% FAST_COVERAGE_SIMULATOR_FUNCTION_V3
%
% Extends v2 with a single additional optimisation:
%   sat_pos_ecef and ue_pos_ecef are cast to single precision immediately after
%   orbit propagation, before the inner UE loop.
%
% Why this helps:
%   double: sat_pos_ecef = 3 × 56 × 2146 × 8 bytes = 2.88 MB
%   single: sat_pos_ecef = 3 × 56 × 2146 × 4 bytes = 1.44 MB
%
%   2.88 MB > 2 MB P-core L2 → every UE iteration must re-fetch from L3.
%   1.44 MB < 2 MB P-core L2 → satellite data stays in L2 after the first UE;
%   subsequent iterations are L2 hits (~8 cycles) rather than L3 hits (~40 cycles).
%
%   L3 thrash onset also doubles: from ~12 simultaneous workers to ~24,
%   so the per-worker time stays low even with more parallel workers.
%
%   Accuracy impact: single precision gives ~7 significant digits.
%   At 1000 km range that is ±0.1 m position error — irrelevant for coverage.
%
% All other logic is identical to v2.
%
% cancel_queue: optional parallel.pool.PollableDataQueue (same as coverage_simulator_function).

    if nargin < 5, cancel_queue = []; end

    fprintf('\n [v3] Fast Simulation: %d Sats, %.1f deg Inc\n', Cfg.Total_sats, Cfg.Inclination);
    tic

    %% 1. Coverage reference latitude
    if isfield(Cfg, 'Lat_range_deg')
        min_lat_cov = min(Cfg.Lat_range_deg);
    else
        min_lat_cov = 0;
    end

    %% 2. Orbit propagation — identical to fast_coverage_simulator_function
    persistent cached_sc last_Cfg
    if nargin < 2 || isempty(reset_cache), reset_cache = false; end
    if reset_cache
        cached_sc = [];  last_Cfg = [];
        fprintf(' [!] Cache reset.\n');
    end

    if use_SGP
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
        r_earth = 6378.14e3;
        if Cfg.WalkerStar
            sats = asymmetrical_walker_star_generation(sc, Cfg.Orbit_height, Cfg.Inclination, ...
                Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Min_elevation_UE, "two-body-keplerian", min_lat_cov);
        else
            sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, Cfg.Inclination, ...
                Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
                Name="S4D", OrbitPropagator="two-body-keplerian");
        end
        [sat_pos_raw, ~, simTimes] = states(sats, "CoordinateFrame", "ECEF");
        sat_pos_ecef = permute(sat_pos_raw, [1, 3, 2]);
    else
        total_duration_sec = seconds(Cfg.StopTime - Cfg.StartTime);
        time_steps_sec     = 0 : Cfg.SampleTime : total_duration_sec;
        simTimes = Cfg.StartTime + seconds(time_steps_sec);
        simTimes.TimeZone = 'UTC';
        if Cfg.WalkerStar
            sat_pos_ecef = fast_walker_star_ecef(Cfg.Orbit_height, Cfg.Inclination, ...
                Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Min_elevation_UE, min_lat_cov, ...
                time_steps_sec, Cfg.StartTime);
        else
            sat_pos_ecef = fast_walker_ecef(Cfg.Orbit_height, Cfg.Inclination, ...
                Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Phasing, time_steps_sec, Cfg.StartTime);
        end
    end

    num_sats = size(sat_pos_ecef, 2);
    nT       = size(sat_pos_ecef, 3);

    %% 2b. Cast to single precision so sat_pos_ecef (1.44 MB) fits in P-core L2 (2 MB)
    %      Each UE iteration reads sat_pos_ecef once; with double it must come from
    %      L3 every time (2.88 MB > L2).  With single it stays in L2 after the first
    %      iteration, cutting ~2146 L3 fetches per call down to ~1.
    sat_pos_ecef = single(sat_pos_ecef);   % 2.88 MB → 1.44 MB

    %% 3. UE positions
    if ~isfield(Cfg, 'Flat_UE_array') || ~isfield(Cfg.Flat_UE_array, 'Lats')
        error('Cfg.Flat_UE_array with Lats/Lons required.');
    end
    UE_lats     = Cfg.Flat_UE_array.Lats(:);
    UE_lons     = Cfg.Flat_UE_array.Lons(:);
    ue_pos_ecef = single(lla2ecef([UE_lats, UE_lons, zeros(numel(UE_lats), 1)]));
    numUEs      = numel(UE_lats);
    Cfg.NumUEs  = numUEs;
    min_el      = Cfg.Min_elevation_UE;

    %% 4. KEY CHANGE: pre-allocate plain matrices, not a struct array
    %
    %  Old code:  UEs(idx).SimData.Range(has_service) = best_ranges(has_service);
    %             → triggers CoW heap alloc + lock on every field write
    %
    %  New code:  best_rng(idx, has_service) = best_ranges(has_service);
    %             → in-place write to a contiguous array, zero extra alloc
    %
    %  For calc_link=false (gridsearch Stage 1/2/3) we only need num_visible.
    %  For calc_link=true  we also need el/az/range/satid to feed the link budget.

    num_vis_mat = zeros(numUEs, nT, 'int16');  % saves 4× memory vs double

    if calc_link
        best_rng = NaN(numUEs, nT, 'single');   % single = 4 bytes vs 8, plenty for km-range
        best_sat = zeros(numUEs, nT, 'int16');
        best_el  = NaN(numUEs, nT, 'single');
        best_az  = NaN(numUEs, nT, 'single');
    end

    %% 5. Inner loop — cancellable serial for-loop (gridsearch compatible)
    %
    %  The parfor path from coverage_simulator_function is NOT used here because:
    %   a) gridsearch always passes cancel_queue → serial path was already forced
    %   b) 32 workers each running their own sim is more efficient than one parfor
    %      over UEs (L3 thrash + heap lock, see analysis in compare_simulator_speed.m)
    %
    %  If you want to run a SINGLE large simulation on the interactive desktop, pass
    %  use_parallel=true via the wrapper at the bottom of this file — it switches the
    %  for to parfor. But gridsearch workers should call this function directly.

    metrics.cancelled          = false;
    metrics.consumed_threshold = Inf;
    consumed_skip              = Inf;

    for idx = 1:numUEs

        % Poll cancel_queue every 100 UEs (identical protocol to coverage_simulator_function)
        if ~isempty(cancel_queue) && mod(idx, 100) == 0
            while true
                [cm, got] = poll(cancel_queue, 0);
                if ~got, break; end
                if isfield(cm, 'skip_above_sats')
                    consumed_skip = min(consumed_skip, cm.skip_above_sats);
                    if Cfg.Total_sats > consumed_skip
                        fprintf('\n[v2 cancel] auto-cancelled at UE %d/%d (%d > thresh %d)\n', ...
                            idx, numUEs, Cfg.Total_sats, consumed_skip);
                        metrics.cancelled          = true;
                        metrics.consumed_threshold = consumed_skip;
                        metrics.worst_coverage_percent = NaN;
                        return;
                    end
                end
                if isfield(cm, 'cancel_detailed') && cm.cancel_detailed
                    fprintf('\n[v2 cancel] Detailed cancelled at UE %d/%d\n', idx, numUEs);
                    metrics.cancelled          = true;
                    metrics.consumed_threshold = consumed_skip;
                    metrics.worst_coverage_percent = NaN;
                    return;
                end
            end
        end

        % -- geometry (same ENU math as before) --
        ue_xyz = ue_pos_ecef(idx, :)';
        dx     = sat_pos_ecef - ue_xyz;           % [3 x num_sats x nT]

        slat = sind(UE_lats(idx));  clat = cosd(UE_lats(idx));
        slon = sind(UE_lons(idx));  clon = cosd(UE_lons(idx));
        R_enu = [-slon,       clon,       0; ...
                 -slat*clon, -slat*slon,  clat; ...
                  clat*clon,  clat*slon,  slat];

        flat    = reshape(dx, 3, []);
        enu_flt = R_enu * flat;
        enu     = reshape(enu_flt, 3, num_sats, nT);

        E = reshape(enu(1,:,:), num_sats, nT);
        N = reshape(enu(2,:,:), num_sats, nT);
        U = reshape(enu(3,:,:), num_sats, nT);

        r_m  = sqrt(E.^2 + N.^2 + U.^2);
        el_m = asind(U ./ r_m);

        valid   = el_m >= min_el;
        nvis    = sum(valid, 1);                         % [1 x nT]

        % -- write to plain matrices (no struct, no CoW) --
        num_vis_mat(idx, :) = int16(nvis);

        if calc_link
            has_srv = nvis > 0;
            if any(has_srv)
                r_tmp = r_m;
                r_tmp(~valid) = Inf;
                [br, bi] = min(r_tmp, [], 1);

                best_rng(idx, has_srv) = single(br(has_srv));

                svc_cols = find(has_srv);
                lin_idxs = bi(has_srv) + (svc_cols - 1) * num_sats;
                best_sat(idx, has_srv) = int16(bi(has_srv));
                best_el(idx,  has_srv) = single(el_m(lin_idxs));

                az_m = atan2d(E, N);
                az_m(az_m < 0) = az_m(az_m < 0) + 360;
                best_az(idx, has_srv) = single(az_m(lin_idxs));
            end
        end

    end  % UE loop

    metrics.consumed_threshold = consumed_skip;
    fprintf('\n [v3] Geometry complete (%.1f sec).\n', toc);

    %% 6. Coverage stats — identical formulas, now from plain matrix
    has_cov      = num_vis_mat >= 1;               % logical [numUEs x nT]
    prob_coverage = 100 * sum(has_cov, 2) / nT;   % [numUEs x 1]

    metrics.worst_coverage_percent = min(prob_coverage);
    metrics.prob_coverage          = prob_coverage;
    metrics.Num_visible            = double(num_vis_mat);
    metrics.minNumberSatellites    = double(min(num_vis_mat, [], 2));
    metrics.meanNumberSatellites   = double(mean(num_vis_mat, 2));
    metrics.throughput_10pct       = NaN;
    metrics.throughput_mean        = NaN;

    if ~calc_link
        %% Build a minimal UEs struct for callers that index into metrics.UEs
        %  (only Lat/Lon/Name + SimData.Num_visible — no heap cost because we
        %   copy already-computed data, not allocating NaN arrays first)
        UEs = struct('Lat', num2cell(UE_lats), ...
                     'Lon', num2cell(UE_lons), ...
                     'Name', strcat('UE', arrayfun(@(k) num2str(k), (1:numUEs)', 'UniformOutput', false)));
        for idx = 1:numUEs
            UEs(idx).SimData = struct( ...
                'Num_visible',   double(num_vis_mat(idx, :)), ...
                'Time',          simTimes, ...
                'SatID',         NaN(1, nT), ...
                'Range',         NaN(1, nT), ...
                'Elevation_deg', NaN(1, nT), ...
                'Azimuth_deg',   NaN(1, nT));
        end
        metrics.UEs = UEs;
        metrics.Cfg = Cfg;
        return;
    end

    %% 7. Full SimData struct (only reached when calc_link=true)
    %     Build once from plain matrices — one alloc per field, no CoW loop.
    UEs(numUEs).Lat = [];
    for idx = 1:numUEs
        UEs(idx).Lat  = UE_lats(idx);
        UEs(idx).Lon  = UE_lons(idx);
        UEs(idx).Name = sprintf('UE%d', idx);
        UEs(idx).SimData = struct( ...
            'Num_visible',   double(num_vis_mat(idx, :)), ...
            'Time',          simTimes, ...
            'SatID',         double(best_sat(idx, :)), ...
            'Range',         double(best_rng(idx, :)), ...
            'Elevation_deg', double(best_el(idx, :)), ...
            'Azimuth_deg',   double(best_az(idx, :)));
    end

    %% 8. Link budget (same batched parfor path as coverage_simulator_function)
    if isfield(Cfg, 'DL')
        tic
        num_workers = Inf;
        batch_size  = 100;
        num_batches = ceil(numUEs / batch_size);

        el_mat_lnk    = double(best_el);
        az_mat_lnk    = double(best_az);
        range_mat_lnk = double(best_rng);
        lat_vec       = UE_lats;
        lon_vec       = UE_lons;

        batch_results = cell(num_batches, 1);
        fprintf('\n [v3] Link budget: %d batches...\n', num_batches);
        parfor (b = 1:num_batches, num_workers)
            i0 = (b-1)*batch_size + 1;
            i1 = min(b*batch_size, numUEs);
            batch_results{b} = link_calc_matrix( ...
                el_mat_lnk(i0:i1,:), az_mat_lnk(i0:i1,:), range_mat_lnk(i0:i1,:), ...
                lat_vec(i0:i1), lon_vec(i0:i1), Cfg.DL, Cfg);
        end

        for b = 1:num_batches
            i0 = (b-1)*batch_size + 1;
            i1 = min(b*batch_size, numUEs);
            ch = batch_results{b};
            for idx = i0:i1
                li = idx - i0 + 1;
                UEs(idx).DL.FSPL       = ch.FSPL(li,:);
                UEs(idx).DL.Total_loss = ch.Total_loss(li,:);
                UEs(idx).DL.SNR        = ch.SNR(li,:);
                UEs(idx).DL.SIR        = ch.SIR(li,:);
                UEs(idx).DL.SINR       = ch.SINR(li,:);
                UEs(idx).DL.Absorption.At = ch.Absorption_At(li,:);
                UEs(idx).DL.PFD_W_MHz  = ch.PFD_W_MHz(li,:);
                UEs(idx).DL.Carrier_density_dBmHz   = ch.Carrier_density_dBmHz(li,:);
                UEs(idx).DL.Noise_density_dBmHz      = ch.Noise_density_dBmHz(li,:);
                UEs(idx).DL.Interference_density_dBmHz = ch.Interference_density_dBmHz(li,:);
                UEs(idx).DL.serving_beam_idx         = ch.serving_beam_idx(li,:);
                UEs(idx).DL.serving_beam_signal_lin  = ch.serving_beam_signal_lin(li,:);
                UEs(idx).DL.interference_lin         = ch.interference_lin(li,:);
                UEs(idx).DL.Frequency  = Cfg.DL.f;
                UEs(idx).DL.Bandwidth  = Cfg.DL.B;
            end
        end

        DL_structs  = [UEs.DL];
        DL_SINR     = vertcat(DL_structs.SINR);
        DL_beam_idx = vertcat(DL_structs.serving_beam_idx);
        SD_structs  = [UEs.SimData];
        SD_SatID    = vertcat(SD_structs.SatID);
        DL_Thpt     = calculate_throughput_matrix(DL_SINR, DL_beam_idx, SD_SatID, Cfg.DL.BeamGrid, Cfg.DL.B, Cfg.Share_bandwidth, Cfg.Modified_shannon);
        for idx = 1:numUEs
            UEs(idx).DL.Throughput = DL_Thpt(idx,:);
        end
        metrics.throughput_10pct = prctile(DL_Thpt(~isnan(DL_Thpt)), 10);
        metrics.throughput_mean  = mean(DL_Thpt(~isnan(DL_Thpt)));
        fprintf('\n [v3] Link budget complete (%.1f sec).\n', toc);
    end

    metrics.UEs                 = UEs;
    metrics.Cfg                 = Cfg;
    SimDataArray                = [UEs.SimData];
    metrics.SimData             = SimDataArray;

end
