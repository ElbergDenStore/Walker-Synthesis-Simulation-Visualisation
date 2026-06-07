function metrics = Constellation_simulator(Cfg, use_parallel, calc_link, use_toolbox, should_cancel)
% CONSTELLATION_SIMULATOR  Satellite coverage + optional link budget.
%
%   metrics = Constellation_simulator(Cfg)
%   metrics = Constellation_simulator(Cfg, use_parallel, calc_link, use_toolbox, should_cancel)
%
% Pipeline
%   propagate orbits -> compute per-UE visibility / best satellite
%                    -> coverage stats -> (optional) link budget
%
% Arguments  (all but Cfg are optional, shown with defaults)
%   Cfg            Struct with StartTime, StopTime, SampleTime, Orbit_height,
%                  Inclination, Num_planes, Sats_per_plane, Total_sats, Phasing,
%                  WalkerStar, Min_elevation_UE, Lat_range_deg,
%                  Flat_UE_array.{Lats,Lons}, and (if calc_link) DL/UL.
%   use_parallel   false. parfor over UEs. Ignored when should_cancel is set.
%   calc_link      false. Compute link budget (much slower).
%   use_toolbox    true.  Satellite Toolbox propagator (cached scenario).
%                         false -> pure-math fast_walker_ecef.
%   should_cancel  @() false. Called every 100 UEs in the inner loop. If it
%                  returns true the run aborts and metrics.cancelled = true.
%                  See CancelToken for the typical gridsearch use.
%
% Returns
%   metrics struct with: worst_coverage_percent, prob_coverage, Num_visible,
%   minNumberSatellites, meanNumberSatellites, throughput_10pct,
%   throughput_mean, UEs, Cfg, cancelled.
%   (Per-UE time series live in UEs(i).SimData; flatten with [metrics.UEs.SimData].)

    %% Defaults
    if nargin < 2 || isempty(use_parallel),  use_parallel  = false;       end
    if nargin < 3 || isempty(calc_link),     calc_link     = false;       end
    if nargin < 4 || isempty(use_toolbox),   use_toolbox   = true;        end
    if nargin < 5 || isempty(should_cancel), should_cancel = @() false;   end

    fprintf('\n Starting Simulation: %d Sats, %.1f deg Inc  [%s%s]\n', ...
        Cfg.Total_sats, Cfg.Inclination, ...
        ternary(use_toolbox, 'toolbox', 'fast-math'), ...
        ternary(calc_link,   ', link budget', ''));

    %% Orbit propagation
    tic;
    [sat_pos_ecef, simTimes] = propagate(Cfg, use_toolbox);
    sat_pos_ecef = single(sat_pos_ecef);   % L2-friendly + faster inner loop

    num_sats = size(sat_pos_ecef, 2);
    nT       = size(sat_pos_ecef, 3);

    %% UE positions
    if ~isfield(Cfg, 'Flat_UE_array') || ~isfield(Cfg.Flat_UE_array, 'Lats')
        error('Cfg.Flat_UE_array.Lats/Lons required.');
    end
    UE_lats = Cfg.Flat_UE_array.Lats(:);
    UE_lons = Cfg.Flat_UE_array.Lons(:);
    numUEs  = numel(UE_lats);
    Cfg.NumUEs  = numUEs;
    ue_pos_ecef = single(lla2ecef([UE_lats, UE_lons, zeros(numUEs, 1)]));
    min_el      = Cfg.Min_elevation_UE;

    %% Result matrices (pre-allocate; written via sliced row writes)
    num_vis_mat = zeros(numUEs, nT, 'int16');
    if calc_link
        best_rng_mat = NaN(numUEs, nT, 'single');
        best_sat_mat = zeros(numUEs, nT, 'int16');
        best_el_mat  = NaN(numUEs, nT, 'single');
        best_az_mat  = NaN(numUEs, nT, 'single');
    else
        best_rng_mat = [];  best_sat_mat = [];
        best_el_mat  = [];  best_az_mat  = [];
    end

    %% Inner UE loop
    %  Serial when cancellation is requested (deterministic poll order).
    %  parfor otherwise iff use_parallel.
    metrics.cancelled = false;
    cancellable       = ~is_default_cancel(should_cancel);

    if cancellable
        for idx = 1:numUEs
            if mod(idx, 100) == 0 && should_cancel()
                metrics.cancelled              = true;
                metrics.worst_coverage_percent = NaN;
                metrics.Num_visible            = [];
                metrics.throughput_10pct       = NaN;
                metrics.throughput_mean        = NaN;
                fprintf('\n[cancel] aborted at UE %d/%d.\n', idx, numUEs);
                return;
            end
            [num_vis_mat(idx,:), r, s, e, a] = compute_ue_row( ...
                ue_pos_ecef(idx,:), sat_pos_ecef, UE_lats(idx), UE_lons(idx), ...
                num_sats, nT, min_el, calc_link);
            if calc_link
                best_rng_mat(idx,:) = r;  best_sat_mat(idx,:) = s;
                best_el_mat(idx,:)  = e;  best_az_mat(idx,:)  = a;
            end
        end
    elseif calc_link
        n_workers = ternary(use_parallel, Inf, 0);
        parfor (idx = 1:numUEs, n_workers)
            [nvis, r, s, e, a] = compute_ue_row( ...
                ue_pos_ecef(idx,:), sat_pos_ecef, UE_lats(idx), UE_lons(idx), ...
                num_sats, nT, min_el, true);
            num_vis_mat(idx,:)  = nvis;
            best_rng_mat(idx,:) = r;  best_sat_mat(idx,:) = s;
            best_el_mat(idx,:)  = e;  best_az_mat(idx,:)  = a;
        end
    else
        n_workers = ternary(use_parallel, Inf, 0);
        parfor (idx = 1:numUEs, n_workers)
            [nvis, ~, ~, ~, ~] = compute_ue_row( ...
                ue_pos_ecef(idx,:), sat_pos_ecef, UE_lats(idx), UE_lons(idx), ...
                num_sats, nT, min_el, false);
            num_vis_mat(idx,:) = nvis;
        end
    end
    fprintf('\nGeometry complete (%.1f sec).\n', toc);

    %% Coverage statistics
    prob_coverage = 100 * sum(num_vis_mat >= 1, 2) ./ nT;
    metrics.worst_coverage_percent = min(prob_coverage);
    metrics.prob_coverage          = prob_coverage;
    metrics.Num_visible            = double(num_vis_mat);
    metrics.minNumberSatellites    = double(min(num_vis_mat, [], 2));
    metrics.meanNumberSatellites   = double(mean(num_vis_mat, 2));
    metrics.throughput_10pct       = NaN;
    metrics.throughput_mean        = NaN;

    %% UEs struct array
    UEs = build_ues_struct(UE_lats, UE_lons, simTimes, num_vis_mat, ...
        best_rng_mat, best_sat_mat, best_el_mat, best_az_mat, calc_link, Cfg, nT);

    %% Optional link budget
    if calc_link && isfield(Cfg, 'DL')
        tic;
        [UEs, metrics.throughput_10pct, metrics.throughput_mean] = ...
            compute_link_budget(UEs, Cfg, use_parallel);
        fprintf('\nLink budget complete (%.1f sec).\n', toc);
    end

    metrics.UEs     = UEs;
    metrics.Cfg     = Cfg;
end

% =========================================================================
%  Propagation
% =========================================================================

function [sat_pos_ecef, simTimes] = propagate(Cfg, use_toolbox)
    if use_toolbox
        [sat_pos_ecef, simTimes] = propagate_toolbox(Cfg);
    else
        [sat_pos_ecef, simTimes] = propagate_fast_math(Cfg);
    end
end

function [sat_pos_ecef, simTimes] = propagate_toolbox(Cfg)
% Cached satelliteScenario + toolbox constellation generators.
    sc   = reuse_or_create_scenario(Cfg);
    sats = generate_constellation(sc, Cfg);
    [raw, ~, simTimes] = states(sats, "CoordinateFrame", "ECEF");
    sat_pos_ecef = permute(raw, [1, 3, 2]);
end

function sc = reuse_or_create_scenario(Cfg)
% Reuse the persistent satelliteScenario container; rebuild only when it
% has been invalidated.  Satellites/ground stations are always rebuilt.
    persistent cached
    if isempty(cached) || ~isvalid(cached)
        cached = satelliteScenario;
    else
        if ~isempty(cached.Satellites),     delete(cached.Satellites);     end
        if ~isempty(cached.GroundStations), delete(cached.GroundStations); end
    end
    sc = cached;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;
end

function sats = generate_constellation(sc, Cfg)
% Walker-Star vs Walker-Delta dispatch.  Hides the asymmetrical-star
% generator's calling convention from the main flow.
    if Cfg.WalkerStar
        min_lat_cov = 0;
        if isfield(Cfg, 'Lat_range_deg'), min_lat_cov = min(Cfg.Lat_range_deg); end
        sats = generate_walker_star_scenario(sc, Cfg.Orbit_height, ...
            Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane, ...
            Cfg.Min_elevation_UE, "two-body-keplerian", min_lat_cov);
    else
        r_earth = 6378.137e3;
        sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, Cfg.Inclination, ...
            Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
            Name="S4D", OrbitPropagator="two-body-keplerian");
    end
end

function [sat_pos_ecef, simTimes] = propagate_fast_math(Cfg)
% Pure-math (no toolbox handle objects).  Used by gridsearch Stage 1 and
% anywhere a Threads pool is needed.
    total_sec      = seconds(Cfg.StopTime - Cfg.StartTime);
    time_steps_sec = 0 : Cfg.SampleTime : total_sec;
    simTimes       = Cfg.StartTime + seconds(time_steps_sec);
    simTimes.TimeZone = 'UTC';
    if Cfg.WalkerStar
        min_lat_cov = 0;
        if isfield(Cfg, 'Lat_range_deg'), min_lat_cov = min(Cfg.Lat_range_deg); end
        sat_pos_ecef = fast_walker_ecef(Cfg.Orbit_height, Cfg.Inclination, ...
            Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Phasing, ...
            time_steps_sec, Cfg.StartTime, true, Cfg.Min_elevation_UE, min_lat_cov);
    else
        sat_pos_ecef = fast_walker_ecef(Cfg.Orbit_height, Cfg.Inclination, ...
            Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Phasing, ...
            time_steps_sec, Cfg.StartTime, false);
    end
end

% =========================================================================
%  Per-UE inner kernel  (vectorised ECEF -> ENU -> elevation)
% =========================================================================

function [nvis_row, rng_row, sat_row, el_row, az_row] = compute_ue_row( ...
    ue_xyz_row, sat_pos_ecef, lat, lon, num_sats, nT, min_el, calc_link)

    ue_xyz = ue_xyz_row(:);                     % 3x1 single
    dx     = sat_pos_ecef - ue_xyz;             % 3 x num_sats x nT

    slat = sind(lat); clat = cosd(lat);
    slon = sind(lon); clon = cosd(lon);
    R_enu = single([ ...
        -slon,        clon,        0;
        -slat*clon,  -slat*slon,   clat;
         clat*clon,   clat*slon,   slat]);

    enu = reshape(R_enu * reshape(dx, 3, []), 3, num_sats, nT);
    E   = reshape(enu(1,:,:), num_sats, nT);
    N   = reshape(enu(2,:,:), num_sats, nT);
    U   = reshape(enu(3,:,:), num_sats, nT);

    r_mat  = sqrt(E.^2 + N.^2 + U.^2);
    el_mat = asind(U ./ r_mat);

    valid    = el_mat >= min_el;
    nvis_row = int16(sum(valid, 1));

    if ~calc_link
        rng_row = [];  sat_row = [];  el_row = [];  az_row = [];
        return;
    end

    rng_row = NaN(1, nT, 'single');
    sat_row = zeros(1, nT, 'int16');
    el_row  = NaN(1, nT, 'single');
    az_row  = NaN(1, nT, 'single');

    has_srv = nvis_row > 0;
    if ~any(has_srv), return; end

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

% =========================================================================
%  UEs struct builder + link budget
% =========================================================================

function UEs = build_ues_struct(UE_lats, UE_lons, simTimes, num_vis_mat, ...
    best_rng_mat, best_sat_mat, best_el_mat, best_az_mat, calc_link, Cfg, nT)

    numUEs = numel(UE_lats);
    UEs(numUEs).Lat = [];   % preallocate

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
            if isfield(Cfg, 'DL'), UEs(idx).DL = empty_link_template(Cfg.DL, nT); end
            if isfield(Cfg, 'UL'), UEs(idx).UL = empty_link_template(Cfg.UL, nT); end
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
    numUEs = numel(UEs);

    SD        = [UEs.SimData];
    el_mat    = single(vertcat(SD.Elevation_deg));
    az_mat    = single(vertcat(SD.Azimuth_deg));
    range_mat = single(vertcat(SD.Range));
    lat_vec   = [UEs.Lat]';
    lon_vec   = [UEs.Lon]';

    n_workers   = ternary(use_parallel, Inf, 0);
    batch_size  = 100;
    num_batches = ceil(numUEs / batch_size);
    batch_results = cell(num_batches, 1);

    fprintf('\nProcessing Link Budget in %d batches...\n', num_batches);
    parfor (b = 1:num_batches, n_workers)
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

    DL_structs     = [UEs.DL];
    DL_SINR        = vertcat(DL_structs.SINR);
    DL_serving_idx = vertcat(DL_structs.serving_beam_idx);
    SimData_SatID  = vertcat(SD.SatID);

    DL_Throughput = calculate_throughput_matrix(DL_SINR, DL_serving_idx, SimData_SatID, ...
        Cfg.DL.BeamGrid, Cfg.DL.B, Cfg.Share_bandwidth, Cfg.Modified_shannon);

    for idx = 1:numUEs
        UEs(idx).DL.Throughput = DL_Throughput(idx, :);
    end

    valid_thpt = DL_Throughput(~isnan(DL_Throughput));
    if isempty(valid_thpt)
        thpt_10pct = NaN;  thpt_mean = NaN;
    else
        thpt_10pct = prctile(valid_thpt, 10);
        thpt_mean  = mean(valid_thpt);
    end
end

% =========================================================================
%  Tiny helpers
% =========================================================================

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function tf = is_default_cancel(fn)
% True if `fn` is the default no-op cancel handle.
    tf = isequal(fn, @() false);
end
