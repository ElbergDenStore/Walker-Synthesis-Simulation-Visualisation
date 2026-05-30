function metrics = constellation_simulator(Cfg, use_parallel, calc_link, cancel_queue, use_toolbox, reset_cache)
% CONSTELLATION_SIMULATOR  Unified satellite coverage + link budget simulator.
%
%   metrics = constellation_simulator(Cfg, use_parallel, calc_link)
%   metrics = constellation_simulator(Cfg, use_parallel, calc_link, cancel_queue)
%   metrics = constellation_simulator(Cfg, use_parallel, calc_link, cancel_queue, use_toolbox)
%   metrics = constellation_simulator(Cfg, use_parallel, calc_link, cancel_queue, use_toolbox, reset_cache)
%
% Replaces the old coverage_simulator_function + fast_coverage_simulator_function pair.
% Combines:
%   * Choice of orbit propagator (Satellite Toolbox SGP/two-body OR pure-math fast_walker_ecef)
%   * Optional link budget computation (calc_link)
%   * Optional cancellation queue for gridsearch (cancel_queue)
%   * Matrix-based inner loop with single-precision satellite positions (v3 optimisation)
%
% Arguments
%   Cfg           Standard config struct (StartTime/StopTime/SampleTime, Total_sats,
%                 Num_planes, Sats_per_plane, Inclination, Phasing, Orbit_height,
%                 Min_elevation_UE, WalkerStar, Lat_range_deg, Flat_UE_array, DL/UL...)
%   use_parallel  true  -> parfor over UEs (only when cancel_queue is empty)
%                 false -> serial for-loop
%   calc_link     true  -> compute full link budget (FSPL, SNR, SINR, throughput, ...)
%                 false -> coverage geometry only (much faster)
%   cancel_queue  Optional parallel.pool.PollableDataQueue. When supplied the inner
%                 loop runs serially and polls every 100 UEs for:
%                   * struct('skip_above_sats', N)  -> auto-cancel if Cfg.Total_sats > N
%                   * struct('cancel_detailed', true) -> explicit cancel
%                 Cancelled runs return early with metrics.cancelled = true.
%   use_toolbox   Default true. true -> Satellite Toolbox; false -> fast_walker_ecef.
%   reset_cache   Default false. Reset the persistent satelliteScenario handle
%                 (only relevant when use_toolbox = true).
%
% Returns metrics struct with fields: worst_coverage_percent, prob_coverage,
% Num_visible, minNumberSatellites, meanNumberSatellites, throughput_10pct,
% throughput_mean, UEs, SimData, Cfg, cancelled, consumed_threshold.

    %% 0. Argument defaults
    if nargin < 4, cancel_queue = []; end
    if nargin < 5 || isempty(use_toolbox),     use_toolbox     = true;  end
    if nargin < 6 || isempty(reset_cache), reset_cache = false; end

    has_cancel_queue = ~isempty(cancel_queue);

    fprintf('\n Starting Simulation: %d Sats, %.1f deg Inc  [propagator: %s%s]\n', ...
        Cfg.Total_sats, Cfg.Inclination, ...
        ternary(use_toolbox, 'toolbox', 'fast-math'), ...
        ternary(calc_link, ', link budget', ''));
    tic

    %% 1. Coverage reference latitude (used by Walker-Star generator only)
    if isfield(Cfg, 'Lat_range_deg')
        min_lat_cov = min(Cfg.Lat_range_deg);
    else
        min_lat_cov = 0;
    end

    %% 2. Orbit propagation
    [sat_pos_ecef, simTimes] = propagate_orbits(Cfg, use_toolbox, reset_cache, min_lat_cov);

    num_sats = size(sat_pos_ecef, 2);
    nT       = size(sat_pos_ecef, 3);

    %% 2b. Cast satellite positions to single precision
    %  Memory: 3 x num_sats x nT x 4 bytes (vs 8 for double).
    %  For a typical Stage-3 grid (56 sats x 2146 steps) this is 1.44 MB instead
    %  of 2.88 MB -> fits in the P-core L2 (2 MB) so the inner UE loop hits L2
    %  instead of L3 on every iteration. Position error ~0.1 m, irrelevant here.
    sat_pos_ecef = single(sat_pos_ecef);

    %% 3. UE positions
    if ~isfield(Cfg, 'Flat_UE_array') || ~isfield(Cfg.Flat_UE_array, 'Lats') || ~isfield(Cfg.Flat_UE_array, 'Lons')
        error('Cfg.Flat_UE_array with Lats/Lons required. Generate UEs before calling constellation_simulator.');
    end
    UE_lats     = Cfg.Flat_UE_array.Lats(:);
    UE_lons     = Cfg.Flat_UE_array.Lons(:);
    numUEs      = numel(UE_lats);
    Cfg.NumUEs  = numUEs;
    ue_pos_ecef = single(lla2ecef([UE_lats, UE_lons, zeros(numUEs, 1)]));

    min_elevation_UE = Cfg.Min_elevation_UE;

    %% 4. Pre-allocate plain result matrices (no struct CoW in the inner loop)
    num_vis_mat = zeros(numUEs, nT, 'int16');
    if calc_link
        best_rng_mat = NaN(numUEs, nT, 'single');
        best_sat_mat = zeros(numUEs, nT, 'int16');
        best_el_mat  = NaN(numUEs, nT, 'single');
        best_az_mat  = NaN(numUEs, nT, 'single');
    else
        % Keep these defined for the post-loop UEs-struct build (some callers
        % rely on metrics.UEs(idx).SimData.Range etc.). Use compact empties so
        % the post-build loop just fills NaNs cheaply.
        best_rng_mat = [];
        best_sat_mat = [];
        best_el_mat  = [];
        best_az_mat  = [];
    end

    %% 5. Inner UE loop
    %  parfor is used only when no cancel_queue is supplied (otherwise the
    %  cancel poll must run serially to be deterministic).
    %
    %  In both paths we accumulate into pre-allocated plain matrices via
    %  sliced row writes (e.g. num_vis_mat(idx,:) = ...). This avoids the
    %  per-iteration struct copy-on-write that dominated the old simulator.
    metrics.cancelled          = false;
    metrics.consumed_threshold = Inf;
    consumed_threshold         = Inf;

    if has_cancel_queue
        % --- SERIAL cancellable loop (gridsearch Stage 3) -------------------
        total_sats_local = Cfg.Total_sats;
        for idx = 1:numUEs
            if mod(idx, 100) == 0
                [cancelled_now, consumed_threshold] = poll_cancel_queue( ...
                    cancel_queue, total_sats_local, consumed_threshold, idx, numUEs);
                if cancelled_now
                    metrics.cancelled              = true;
                    metrics.consumed_threshold     = consumed_threshold;
                    metrics.worst_coverage_percent = NaN;
                    metrics.Num_visible            = [];
                    metrics.SimData                = [];
                    metrics.throughput_10pct       = NaN;
                    metrics.throughput_mean        = NaN;
                    return;
                end
            end

            [nvis_row, rng_row, sat_row, el_row, az_row] = compute_ue_row( ...
                ue_pos_ecef(idx,:), sat_pos_ecef, UE_lats(idx), UE_lons(idx), ...
                num_sats, nT, min_elevation_UE, calc_link);

            num_vis_mat(idx, :) = nvis_row;
            if calc_link
                best_rng_mat(idx, :) = rng_row;
                best_sat_mat(idx, :) = sat_row;
                best_el_mat(idx,  :) = el_row;
                best_az_mat(idx,  :) = az_row;
            end
        end
    else
        % --- Optionally parallel loop (no cancellation) ---------------------
        if use_parallel, num_workers = Inf; else, num_workers = 0; end

        if calc_link
            parfor (idx = 1:numUEs, num_workers)
                [nvis_row, rng_row, sat_row, el_row, az_row] = compute_ue_row( ...
                    ue_pos_ecef(idx,:), sat_pos_ecef, UE_lats(idx), UE_lons(idx), ...
                    num_sats, nT, min_elevation_UE, true);
                num_vis_mat(idx, :)  = nvis_row;
                best_rng_mat(idx, :) = rng_row;
                best_sat_mat(idx, :) = sat_row;
                best_el_mat(idx,  :) = el_row;
                best_az_mat(idx,  :) = az_row;
            end
        else
            parfor (idx = 1:numUEs, num_workers)
                [nvis_row, ~, ~, ~, ~] = compute_ue_row( ...
                    ue_pos_ecef(idx,:), sat_pos_ecef, UE_lats(idx), UE_lons(idx), ...
                    num_sats, nT, min_elevation_UE, false);
                num_vis_mat(idx, :) = nvis_row;
            end
        end
    end

    metrics.consumed_threshold = consumed_threshold;
    fprintf('\nGeometry complete (%.1f sec).\n', toc);

    %% 6. Coverage statistics (from matrices)
    has_cov              = num_vis_mat >= 1;
    prob_coverage        = 100 * sum(has_cov, 2) ./ nT;
    minNumberSatellites  = double(min(num_vis_mat, [], 2));
    meanNumberSatellites = double(mean(num_vis_mat, 2));

    metrics.worst_coverage_percent = min(prob_coverage);
    metrics.prob_coverage          = prob_coverage;
    metrics.Num_visible            = double(num_vis_mat);
    metrics.minNumberSatellites    = minNumberSatellites;
    metrics.meanNumberSatellites   = meanNumberSatellites;
    metrics.throughput_10pct       = NaN;
    metrics.throughput_mean        = NaN;

    %% 7. Build UEs struct array from matrices (one allocation per field)
    UEs = build_ues_struct(UE_lats, UE_lons, simTimes, num_vis_mat, ...
        best_rng_mat, best_sat_mat, best_el_mat, best_az_mat, calc_link, Cfg, nT);

    %% 8. Optional link budget
    if calc_link && isfield(Cfg, 'DL')
        tic
        [UEs, metrics.throughput_10pct, metrics.throughput_mean] = ...
            compute_link_budget(UEs, Cfg, use_parallel);
        fprintf('\nLink budget complete (%.1f sec).\n', toc);
    end

    %% 9. Package
    metrics.UEs     = UEs;
    metrics.Cfg     = Cfg;
    metrics.SimData = [UEs.SimData];
end

% =========================================================================
% Local helpers
% =========================================================================

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function [sat_pos_ecef, simTimes] = propagate_orbits(Cfg, use_toolbox, reset_cache, min_lat_cov)
% Returns [3 x num_sats x nT] double position array + simTimes datetime row.
    persistent cached_sc
    if reset_cache, cached_sc = []; end

    if use_toolbox
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
            sats = asymmetrical_walker_star_generation(sc, Cfg.Orbit_height, ...
                Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane, ...
                Cfg.Min_elevation_UE, "two-body-keplerian", min_lat_cov);
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
                Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Min_elevation_UE, ...
                min_lat_cov, time_steps_sec, Cfg.StartTime);
        else
            sat_pos_ecef = fast_walker_ecef(Cfg.Orbit_height, Cfg.Inclination, ...
                Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Phasing, ...
                time_steps_sec, Cfg.StartTime);
        end
    end
end

function [cancelled, consumed_threshold] = poll_cancel_queue(q, total_sats, consumed_threshold, idx, numUEs)
% Drain pending messages from the cancel queue. Returns cancelled=true if either
%   * a skip_above_sats threshold has been set below the current Total_sats, or
%   * an explicit cancel_detailed message has arrived.
    cancelled = false;
    while true
        [msg, got] = poll(q, 0);
        if ~got, break; end
        if isfield(msg, 'skip_above_sats')
            consumed_threshold = min(consumed_threshold, msg.skip_above_sats);
            if total_sats > consumed_threshold
                fprintf('\n[cancel] auto-cancelled at UE %d/%d (%d sats > threshold %d).\n', ...
                    idx, numUEs, total_sats, consumed_threshold);
                cancelled = true;
                return;
            end
        end
        if isfield(msg, 'cancel_detailed') && msg.cancel_detailed
            fprintf('\n[cancel] explicit cancel at UE %d/%d.\n', idx, numUEs);
            cancelled = true;
            return;
        end
    end
end

function [nvis_row, rng_row, sat_row, el_row, az_row] = compute_ue_row( ...
    ue_xyz_row, sat_pos_ecef, lat, lon, num_sats, nT, min_el, calc_link)
% Compute one UE's row of visibility / best-satellite data.
%   ue_xyz_row  1x3 single (ECEF position of this UE)
%   Returns row vectors (1 x nT) suitable for sliced parfor writes.

    ue_xyz = ue_xyz_row(:);                          % 3x1 single
    dx     = sat_pos_ecef - ue_xyz;                  % [3 x num_sats x nT] single

    slat = sind(lat); clat = cosd(lat);
    slon = sind(lon); clon = cosd(lon);
    R_enu = single([ ...
        -slon,        clon,        0;
        -slat*clon,  -slat*slon,   clat;
         clat*clon,   clat*slon,   slat]);

    flat    = reshape(dx, 3, []);
    enu_flt = R_enu * flat;
    enu     = reshape(enu_flt, 3, num_sats, nT);

    E = reshape(enu(1,:,:), num_sats, nT);
    N = reshape(enu(2,:,:), num_sats, nT);
    U = reshape(enu(3,:,:), num_sats, nT);

    r_mat  = sqrt(E.^2 + N.^2 + U.^2);
    el_mat = asind(U ./ r_mat);

    valid    = el_mat >= min_el;
    nvis_row = int16(sum(valid, 1));                 % 1 x nT

    if calc_link
        rng_row = NaN(1, nT, 'single');
        sat_row = zeros(1, nT, 'int16');
        el_row  = NaN(1, nT, 'single');
        az_row  = NaN(1, nT, 'single');

        has_srv = nvis_row > 0;
        if any(has_srv)
            r_tmp = r_mat;
            r_tmp(~valid) = Inf;
            [br, bi] = min(r_tmp, [], 1);

            rng_row(has_srv) = br(has_srv);
            sat_row(has_srv) = int16(bi(has_srv));

            svc_cols = find(has_srv);
            lin_idxs = bi(has_srv) + (svc_cols - 1) * num_sats;
            el_row(has_srv) = el_mat(lin_idxs);

            az_mat = atan2d(E, N);
            az_mat(az_mat < 0) = az_mat(az_mat < 0) + 360;
            az_row(has_srv) = az_mat(lin_idxs);
        end
    else
        rng_row = [];
        sat_row = [];
        el_row  = [];
        az_row  = [];
    end
end

function UEs = build_ues_struct(UE_lats, UE_lons, simTimes, num_vis_mat, ...
    best_rng_mat, best_sat_mat, best_el_mat, best_az_mat, calc_link, Cfg, nT)
% Build the UEs struct array once from already-computed matrices.
% Avoids the per-iteration CoW that made the old loop slow.

    numUEs = numel(UE_lats);
    UEs(numUEs).Lat = [];   % allocate

    for idx = 1:numUEs
        UEs(idx).Lat  = UE_lats(idx);
        UEs(idx).Lon  = UE_lons(idx);
        UEs(idx).Name = sprintf('UE%d', idx);

        if calc_link
            UEs(idx).SimData = struct( ...
                'Time',          simTimes, ...
                'SatID',         double(best_sat_mat(idx, :)), ...
                'Range',         double(best_rng_mat(idx, :)), ...
                'Elevation_deg', double(best_el_mat(idx, :)), ...
                'Azimuth_deg',   double(best_az_mat(idx, :)), ...
                'Num_visible',   double(num_vis_mat(idx, :)));
            % DL/UL placeholders so compute_link_budget can fill them in.
            if isfield(Cfg, 'DL')
                UEs(idx).DL = empty_link_template(Cfg.DL, nT);
            end
            if isfield(Cfg, 'UL')
                UEs(idx).UL = empty_link_template(Cfg.UL, nT);
            end
        else
            UEs(idx).SimData = struct( ...
                'Time',          simTimes, ...
                'SatID',         NaN(1, nT), ...
                'Range',         NaN(1, nT), ...
                'Elevation_deg', NaN(1, nT), ...
                'Azimuth_deg',   NaN(1, nT), ...
                'Num_visible',   double(num_vis_mat(idx, :)));
        end
    end
end

function L = empty_link_template(linkCfg, nT)
    L = struct( ...
        'Frequency',                   linkCfg.f, ...
        'Bandwidth',                   linkCfg.B, ...
        'FSPL',                        NaN(1, nT), ...
        'Absorption',                  struct('At', NaN(1, nT)), ...
        'Total_loss',                  NaN(1, nT), ...
        'PFD_W_MHz',                   NaN(1, nT), ...
        'Noise_density_dBmHz',         NaN(1, nT), ...
        'Carrier_density_dBmHz',       NaN(1, nT), ...
        'Interference_density_dBmHz',  NaN(1, nT), ...
        'SNR',                         NaN(1, nT), ...
        'SIR',                         NaN(1, nT), ...
        'SINR',                        NaN(1, nT), ...
        'serving_beam_idx',            NaN(1, nT), ...
        'serving_beam_signal_lin',     NaN(1, nT), ...
        'interference_lin',            NaN(1, nT), ...
        'Throughput',                  NaN(1, nT));
end

function [UEs, thpt_10pct, thpt_mean] = compute_link_budget(UEs, Cfg, use_parallel)
% Same batched-parfor link budget path as the original, but driven from the
% already-built UEs struct (no per-loop CoW writes).
    numUEs = numel(UEs);

    SimDataArray = [UEs.SimData];
    el_mat    = single(vertcat(SimDataArray.Elevation_deg));
    az_mat    = single(vertcat(SimDataArray.Azimuth_deg));
    range_mat = single(vertcat(SimDataArray.Range));
    lat_vec   = [UEs.Lat]';
    lon_vec   = [UEs.Lon]';

    if use_parallel, num_workers = Inf; else, num_workers = 0; end

    batch_size  = 100;
    num_batches = ceil(numUEs / batch_size);
    batch_results = cell(num_batches, 1);

    fprintf('\nProcessing Link Budget in %d batches...\n', num_batches);
    parfor (b = 1:num_batches, num_workers)
        i0 = (b - 1) * batch_size + 1;
        i1 = min(b * batch_size, numUEs);
        batch_results{b} = link_calc_matrix( ...
            el_mat(i0:i1, :), az_mat(i0:i1, :), range_mat(i0:i1, :), ...
            lat_vec(i0:i1), lon_vec(i0:i1), Cfg.DL, Cfg);
    end

    for b = 1:num_batches
        i0 = (b - 1) * batch_size + 1;
        i1 = min(b * batch_size, numUEs);
        ch = batch_results{b};
        li = 0;
        for idx = i0:i1
            li = li + 1;
            UEs(idx).DL.FSPL                       = ch.FSPL(li, :);
            UEs(idx).DL.Total_loss                 = ch.Total_loss(li, :);
            UEs(idx).DL.PFD_W_MHz                  = ch.PFD_W_MHz(li, :);
            UEs(idx).DL.Noise_density_dBmHz        = ch.Noise_density_dBmHz(li, :);
            UEs(idx).DL.Carrier_density_dBmHz      = ch.Carrier_density_dBmHz(li, :);
            UEs(idx).DL.Interference_density_dBmHz = ch.Interference_density_dBmHz(li, :);
            UEs(idx).DL.SNR                        = ch.SNR(li, :);
            UEs(idx).DL.SIR                        = ch.SIR(li, :);
            UEs(idx).DL.SINR                       = ch.SINR(li, :);
            UEs(idx).DL.serving_beam_idx           = ch.serving_beam_idx(li, :);
            UEs(idx).DL.serving_beam_signal_lin    = ch.serving_beam_signal_lin(li, :);
            UEs(idx).DL.interference_lin           = ch.interference_lin(li, :);
            UEs(idx).DL.Absorption.At              = ch.Absorption_At(li, :);
        end
    end

    DL_structs        = [UEs.DL];
    DL_SINR           = vertcat(DL_structs.SINR);
    DL_serving_idx    = vertcat(DL_structs.serving_beam_idx);
    SD_structs        = [UEs.SimData];
    SimData_SatID     = vertcat(SD_structs.SatID);

    DL_Throughput = calculate_throughput_matrix(DL_SINR, DL_serving_idx, SimData_SatID, ...
        Cfg.DL.BeamGrid, Cfg.DL.B, Cfg.Share_bandwidth, Cfg.Modified_shannon);

    for idx = 1:numUEs
        UEs(idx).DL.Throughput = DL_Throughput(idx, :);
    end

    valid_thpt = DL_Throughput(~isnan(DL_Throughput));
    if isempty(valid_thpt)
        thpt_10pct = NaN;
        thpt_mean  = NaN;
    else
        thpt_10pct = prctile(valid_thpt, 10);
        thpt_mean  = mean(valid_thpt);
    end
end
