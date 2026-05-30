function [best_params, all_candidates, out_dir] = gridsearch(master_config, orbit_height_km, min_sats, base_out_dir)
    if nargin < 4 || isempty(base_out_dir)
        base_out_dir = fullfile('simulation_output', 'gridsearch_runs');
    end
    %% 1. Build the Ascending Grid
    % Walker Star mode: phasing fixed at p/2 (can be non-integer); no phasing loop.
    % Walker Delta mode: phase factor f = 0:(p-1) gives p values per plane count.
    walker_star_mode = isfield(master_config, 'WalkerStar') && master_config.WalkerStar;
    if walker_star_mode
        Num_constellations = numel(master_config.Num_Planes) * numel(master_config.Sats_Plane) * numel(master_config.Inc_vec);
    else
        Num_constellations = sum(master_config.Num_Planes) * numel(master_config.Sats_Plane) * numel(master_config.Inc_vec);
    end
    grid_data = zeros(Num_constellations,5);
    i = 1;
    for p = master_config.Num_Planes
        for s = master_config.Sats_Plane
            if walker_star_mode
                phasing_vec = p/2;
            else
                phasing_vec = 0:(p-1);
            end
            for f = phasing_vec
                for inc = master_config.Inc_vec
                    grid_data(i,:) = [p, s, inc, f, p*s];
                    i = i + 1;
                end
            end
        end
    end

    search_grid = array2table(grid_data, 'VariableNames', ...
        {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing', 'Total_sats'});

    %% 2. Smart Filters
    isValidTarget = search_grid.Total_sats >= min_sats;
    search_grid = search_grid(isValidTarget, :);

    % Sort from cheapest to most expensive
    search_grid = sortrows(search_grid, 'Total_sats', 'ascend');

    fprintf('\n=== Starting Ascending Grid Search ===\n');
    fprintf('Testing %d valid architectures from %d satellites...\n\n', height(search_grid), min_sats);

    %% 3. Setup Parallel Environment and Queue
    best_params = [];
    all_candidates = table();
    target_num_candidates = master_config.Target_num_candidates;

    pool = gcp('nocreate');
    if isempty(pool), pool = parpool(); end
    num_workers = pool.NumWorkers;
    total_runs = height(search_grid);

    fprintf('Using %d parallel workers for batch processing...\n', num_workers);

    % --- DASHBOARD & PROFILING TRACKERS ---
    t_faster_all = NaN(total_runs, 1);
    t_fast_all = NaN(total_runs, 1);
    t_detailed_all = NaN(total_runs, 1);
    worker_run_counts = zeros(num_workers, 1);
    worker_total_math_time = zeros(num_workers, 1);

    w_states = repmat("Idle", 1, num_workers);
    w_runs = zeros(1, num_workers);
    w_last_msg_s = zeros(1, num_workers);
    event_log = repmat({''}, 1, 10);

    if isfield(master_config, 'Worker_stall_timeout_s')
        stall_timeout_s = master_config.Worker_stall_timeout_s;
    else
        stall_timeout_s = 3600; % 1 hour — low-altitude runs with many sats can take >500s
    end
    % Grace period after cancel_detailed is sent before the worker is force-killed.
    % constellation_simulator polls every 100 UEs; 5 min is ample.
    cancel_grace_s  = 300;
    cancel_sent_at  = zeros(1, num_workers); % wall-time when cancel was sent (0 = not sent)

    evaluated_coverage = NaN(total_runs, 1);
    detailed_coverage = NaN(total_runs, 1);
    status_flags = zeros(total_runs, 1); % 0 = Invalid, 1 = Candidate, 2 = Best
    is_completed = false(total_runs, 1);

    min_sats_found = Inf;
    candidates_found = 0;
    runs_completed = 0;
    worst_accepted_sats = NaN;

    q = parallel.pool.PollableDataQueue;
    futures = parallel.FevalFuture.empty(1, 0);
    worker_assigned_indices = cell(1, num_workers); % tracks which run indices each worker owns
    worker_input_queues = cell(1, num_workers);     % back-channel queue (client -> worker)
    runs_skipped = 0;                               % counts runs dropped due to worker death
    runs_skipped_threshold = 0;                     % runs skipped because they exceed worst-accepted sats
    threshold_broadcast = false;                    % becomes true once skip-threshold is sent
    skip_above_sats = Inf;                          % cached threshold for late-registering workers

    %% 4. Submit Interleaved Chunks (Glass Cockpit model)
    fprintf('\nSubmitting interleaved chunks to %d workers...\n\n', num_workers);

    for w = 1:num_workers
        % Interleave indices so all workers start from cheaper architectures.
        indices = w:num_workers:total_runs;
        if ~isempty(indices)
            worker_assigned_indices{w} = indices;
            worker_grid = search_grid(indices, :);
            futures(w) = parfeval(pool, @evaluate_chunk, 0, ...
                q, master_config, orbit_height_km, worker_grid, indices, w);
            w_states(w) = "Queued";
            w_runs(w) = indices(1);
        end
    end

    %% 5. Synchronous Client Loop (The Glass Cockpit)
    t_run_start = tic;
    t_last_update = tic;

    while true
        [data, gotMsg] = poll(q, 0.01);

        if gotMsg
            wid = data.worker_id;
            w_last_msg_s(wid) = toc(t_run_start);

            if isfield(data, 'is_registration') && data.is_registration
                worker_input_queues{wid} = data.queue_handle;
                % If threshold was already broadcast before this worker registered, send it now.
                if threshold_broadcast
                    try
                        send(worker_input_queues{wid}, struct('skip_above_sats', skip_above_sats));
                    catch
                        worker_input_queues{wid} = [];
                    end
                end
                continue;
            end

            w_runs(wid) = data.run_idx;

            if isfield(data, 'is_worker_error') && data.is_worker_error
                error('\n[!] WORKER %d ERROR ON RUN %d: %s', wid, data.run_idx, data.error_message);
            end

            if isfield(data, 'is_skipped') && data.is_skipped
                runs_completed = runs_completed + 1;
                runs_skipped_threshold = runs_skipped_threshold + 1;
                is_completed(data.run_idx) = true;
                w_states(wid) = "Skip";
                continue;
            end

            if isfield(data, 'is_heartbeat') && data.is_heartbeat
                w_states(wid) = data.state;

                if strcmp(data.state, "DETAILED")
                    new_event = sprintf('[W%02d] STARTING DETAILED: Run %d | %s', ...
                        wid, data.run_idx, data.arch_msg);
                    event_log = [event_log(2:end), {new_event}];
                end
            else
                runs_completed = runs_completed + 1;
                idx = data.run_idx;
                is_completed(idx) = true;

                t_faster_all(idx) = data.t_faster;
                if data.t_fast > 0, t_fast_all(idx) = data.t_fast; end
                if data.t_detailed > 0, t_detailed_all(idx) = data.t_detailed; end

                worker_run_counts(wid) = worker_run_counts(wid) + 1;
                worker_total_math_time(wid) = worker_total_math_time(wid) + data.t_total;

                w_states(wid) = "Ultra";
                evaluated_coverage(idx) = data.faster_cov;
                if isfield(data, 'detailed_cov') && ~isnan(data.detailed_cov)
                    detailed_coverage(idx) = data.detailed_cov;
                end

                if data.is_candidate
                    new_event = sprintf('[W%02d] Run %d: *** %s ***', wid, idx, data.msg);
                    event_log = [event_log(2:end), {new_event}];

                    candidates_found = candidates_found + 1;
                    status_flags(idx) = 1;
                    temp_params = search_grid(idx, :);
                    temp_params.Orbit_height_km = orbit_height_km;
                    all_candidates = [all_candidates; temp_params];

                    if temp_params.Total_sats < min_sats_found
                        min_sats_found = temp_params.Total_sats;
                        best_params = temp_params;
                    end
                end

                % --- DETERMINISTIC STOPPING / THRESHOLD BROADCAST ---
                if candidates_found >= target_num_candidates
                    sorted_candidates = sortrows(all_candidates, 'Total_sats', 'ascend');
                    worst_accepted_sats = sorted_candidates.Total_sats(target_num_candidates);

                    first_pending_idx = find(~is_completed, 1, 'first');
                    if isempty(first_pending_idx)
                        break;
                    end

                    % Broadcast (or update) the skip threshold to all workers so they
                    % skip remaining runs above the worst-accepted satellite count.
                    % Workers that pass the threshold continue normally; if everything
                    % left in a worker's chunk is above the threshold, it exits early.
                    if ~threshold_broadcast || worst_accepted_sats < skip_above_sats
                        skip_above_sats = worst_accepted_sats;
                        for wq = 1:numel(worker_input_queues)
                            if ~isempty(worker_input_queues{wq})
                                try
                                    send(worker_input_queues{wq}, struct('skip_above_sats', skip_above_sats));
                                catch
                                    % Worker died before receiving the threshold — clean up its handle
                                    % so we don't try again; the future-state check will mark it dead.
                                    worker_input_queues{wq} = [];
                                end
                            end
                        end
                        if ~threshold_broadcast
                            fprintf('\n--- Skip threshold broadcast: drop runs with > %d sats (workers idle once chunk drained) ---\n', skip_above_sats);
                        end
                        threshold_broadcast = true;
                    end
                end
            end
        end

        % Refresh dashboard every 0.5 s.
        if toc(t_last_update) > 0.5
            print_dashboard(w_states, w_runs, runs_completed, total_runs, toc(t_run_start), event_log, candidates_found, target_num_candidates);
            t_last_update = tic;
        end

        f_states = {futures.State};
        terminal_states = {'finished', 'unavailable', 'failed', 'cancelled'};
        if all(ismember(f_states, terminal_states))
            % Drain any remaining queue messages (e.g. error reports sent just before worker died)
            while true
                [drain_data, got_drain] = poll(q, 0.001);
                if ~got_drain, break; end
                if isfield(drain_data, 'is_worker_error') && drain_data.is_worker_error
                    error('\n[!] WORKER %d ERROR ON RUN %d: %s', drain_data.worker_id, drain_data.run_idx, drain_data.error_message);
                end
            end
            % Check for crashes — skip dead workers' remaining runs instead of aborting
            for fe = 1:numel(futures)
                is_dead = ~strcmp(futures(fe).State, 'finished');
                if is_dead && ~isempty(worker_assigned_indices{fe})
                    msg = futures(fe).State;
                    if ~isempty(futures(fe).Error), msg = futures(fe).Error.message; end
                    incomplete_idx = worker_assigned_indices{fe}(~is_completed(worker_assigned_indices{fe}));
                    is_completed(incomplete_idx) = true;
                    runs_skipped = runs_skipped + numel(incomplete_idx);
                    worker_assigned_indices{fe} = []; % prevent double-processing
                    fprintf('\n[!] WORKER %d DIED: %s. Skipped %d remaining runs.\n', fe, msg, numel(incomplete_idx));
                end
            end
            break;
        end

        if any(~ismember(f_states, {'running', 'queued', 'finished'}))
            for e = 1:numel(futures)
                is_dead = ~ismember(f_states{e}, {'running', 'queued', 'finished'});
                if is_dead && ~isempty(worker_assigned_indices{e})
                    msg = f_states{e};
                    if ~isempty(futures(e).Error), msg = futures(e).Error.message; end
                    incomplete_idx = worker_assigned_indices{e}(~is_completed(worker_assigned_indices{e}));
                    is_completed(incomplete_idx) = true;
                    runs_skipped = runs_skipped + numel(incomplete_idx);
                    worker_assigned_indices{e} = []; % prevent double-processing
                    fprintf('\n[!] WORKER %d DIED (%s). Skipped %d remaining runs.\n', e, msg, numel(incomplete_idx));
                end
            end
        end

        % Cooperatively cancel DETAILED workers whose constellation is above the
        % skip threshold — sends a message via q_in; the worker polls it every
        % 50 UEs inside constellation_simulator and returns early.
        if threshold_broadcast
            for w = 1:length(futures)
                if strcmp(w_states(w), "DETAILED") && ~isempty(worker_input_queues{w})
                    cur_idx = w_runs(w);
                    if cur_idx >= 1 && cur_idx <= total_runs && ...
                            search_grid.Total_sats(cur_idx) > skip_above_sats
                        fprintf('\n[i] W%02d cancel requested: %d sats > threshold %d\n', ...
                            w, search_grid.Total_sats(cur_idx), skip_above_sats);
                        try
                            send(worker_input_queues{w}, struct('cancel_detailed', true));
                        catch
                            worker_input_queues{w} = [];
                        end
                        w_states(w) = "Cancelling";
                        cancel_sent_at(w) = toc(t_run_start);
                    end
                end
            end
        end

        for w = 1:length(futures)
            if strcmp(f_states{w}, 'running')
                % Force-kill workers that ignored a cancel signal for too long.
                if strcmp(w_states(w), "Cancelling") && cancel_sent_at(w) > 0
                    if (toc(t_run_start) - cancel_sent_at(w)) > cancel_grace_s
                        cancel(futures(w));
                        incomplete_idx = worker_assigned_indices{w}(~is_completed(worker_assigned_indices{w}));
                        is_completed(incomplete_idx) = true;
                        runs_skipped = runs_skipped + numel(incomplete_idx);
                        worker_assigned_indices{w} = [];
                        worker_input_queues{w} = [];
                        cancel_sent_at(w) = 0;
                        fprintf('\n[!] W%02d ignored cancel for >%.0fs — force-killed. Skipped %d runs.\n', ...
                            w, cancel_grace_s, numel(incomplete_idx));
                        w_states(w) = "Killed";
                    end
                % Regular stall detection for non-cancelling workers.
                elseif (toc(t_run_start) - w_last_msg_s(w)) > stall_timeout_s
                    error(['[FATAL] Worker %d stalled for >%.0fs (state=%s, run=%d). ' ...
                        'Crashing to preserve optimality guarantee. ' ...
                        'Increase Worker_stall_timeout_s if runs legitimately take this long.'], ...
                        w, stall_timeout_s, w_states(w), w_runs(w));
                end
            end
        end
    end

    %% 6. End-of-loop cleanup
    cancel(futures);
    run_time = toc(t_run_start);

    generate_profiling_report(run_time, runs_completed, num_workers, ...
        t_faster_all, t_fast_all, t_detailed_all, worker_run_counts, worker_total_math_time);
    if runs_skipped > 0
        fprintf('[!] %d runs skipped due to worker failure(s) (%.1f%% of filtered grid).\n', ...
            runs_skipped, 100*runs_skipped/max(total_runs,1));
    end
    if runs_skipped_threshold > 0
        fprintf('[i] %d runs skipped because total_sats exceeded worst-accepted (%d).\n', ...
            runs_skipped_threshold, skip_above_sats);
    end

    out_dir = ''; % initialised here; assigned in section 7 when a candidate is found
    if isempty(best_params)
        fprintf('\n[!] GRID SEARCH EXHAUSTED [!]\n');
        fprintf('No constellation achieved 99.999%% coverage within limits.\n');
        return;
    end

    best_sats_count = min(search_grid.Total_sats(status_flags > 0));
    if ~isempty(best_sats_count)
        best_indices = find(search_grid.Total_sats == best_sats_count & status_flags > 0);
        status_flags(best_indices) = 2;
    end

    %% 7. Create Output Directory (report + plots)
    date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    folder_name = sprintf('%.0f_%d_%s', orbit_height_km, best_params.Total_sats, date_str);
    out_dir = fullfile(base_out_dir, folder_name);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    report_path = fullfile(out_dir, 'deep_profile_report.txt');
    write_profiling_report_to_file(report_path, master_config, run_time, runs_completed, num_workers, ...
        t_faster_all, t_fast_all, t_detailed_all, worker_run_counts, worker_total_math_time, ...
        best_params, all_candidates, search_grid, status_flags, Num_constellations, total_runs);

    %% 8. Save Plot Data
    t1_valid = t_faster_all(~isnan(t_faster_all));
    t2_valid = t_fast_all(~isnan(t_fast_all));
    t3_valid = t_detailed_all(~isnan(t_detailed_all));
    n1 = numel(t1_valid); n2 = numel(t2_valid); n3 = numel(t3_valid);

    plot_data.search_grid           = search_grid;
    plot_data.status_flags          = status_flags;
    plot_data.evaluated_coverage    = evaluated_coverage;
    plot_data.detailed_coverage     = detailed_coverage;
    plot_data.orbit_height_km       = orbit_height_km;
    plot_data.worst_accepted_sats   = worst_accepted_sats;
    plot_data.all_candidates        = all_candidates;
    plot_data.target_num_candidates = target_num_candidates;
    % Multi-stage profiling stats (consumed by summarize_sweep)
    plot_data.n_full_grid           = Num_constellations;
    plot_data.n_total_grid          = total_runs;
    plot_data.n_stage1              = n1;
    plot_data.n_stage2              = n2;
    plot_data.n_stage3              = n3;
    plot_data.avg_t1_s              = sum(t1_valid) / max(n1, 1);
    plot_data.avg_t2_s              = sum(t2_valid) / max(n2, 1);
    plot_data.avg_t3_s              = sum(t3_valid) / max(n3, 1);
    plot_data.wall_time_s           = run_time;
    plot_data.num_workers           = num_workers;
    plot_data.candidates_found      = candidates_found;
    plot_data.runs_skipped          = runs_skipped;
    save(fullfile(out_dir, 'plot_data.mat'), 'plot_data');
    fprintf('Plot data saved. Regenerate plots with: plot_gridsearch(''%s'')\n', out_dir);
end

% =========================================================================
% WORKER & HELPER FUNCTIONS
% =========================================================================
function evaluate_chunk(q, master_config, orbit_height_km, worker_grid, original_indices, worker_id)
    % Force workers to one thread to avoid oversubscription.
    maxNumCompThreads(1);

    % Back-channel queue so the client can broadcast a skip threshold.
    % worker_grid is a sorted-ascending stride of search_grid, so once a row
    % exceeds the threshold every remaining row does too -> safe to early-exit.
    q_in = parallel.pool.PollableDataQueue;
    send(q, struct('worker_id', worker_id, 'is_registration', true, 'queue_handle', q_in));
    skip_above_sats = Inf;

    for i = 1:height(worker_grid)
        % Drain any pending threshold updates (keep the most recent / smallest).
        while true
            [msg, gotMsg] = poll(q_in, 0);
            if ~gotMsg, break; end
            if isfield(msg, 'skip_above_sats')
                skip_above_sats = min(skip_above_sats, msg.skip_above_sats);
            end
        end

        if worker_grid.Total_sats(i) > skip_above_sats
            % All remaining rows are >= this one (sorted ascending) -> skip them all.
            for j = i:height(worker_grid)
                send(q, struct('worker_id', worker_id, ...
                    'run_idx', original_indices(j), ...
                    'is_skipped', true));
            end
            return;
        end

        local_config.Orbit_height_m = orbit_height_km * 1e3;
        local_config.Num_planes = worker_grid.Num_planes(i);
        local_config.Sats_per_plane = worker_grid.Sats_per_plane(i);
        local_config.Inclination = worker_grid.Inclination(i);
        local_config.Phasing = worker_grid.Phasing(i);
        local_config.Total_sats = worker_grid.Total_sats(i);

        try
            res = run_single_evaluation(q, local_config, master_config, original_indices(i), worker_id, q_in, skip_above_sats);
            % Pull back any threshold the worker consumed from q_in during the run
            if isfield(res, 'updated_skip_above_sats')
                skip_above_sats = min(skip_above_sats, res.updated_skip_above_sats);
            end
            res.is_heartbeat = false;
            send(q, res);
        catch ME
            send(q, struct('worker_id', worker_id, ...
                'run_idx', original_indices(i), ...
                'is_worker_error', true, ...
                'error_message', ME.message));
            rethrow(ME);
        end
    end
end

function result = run_single_evaluation(q, local_config, master_config, run_idx, worker_id, q_in, skip_above_sats)
    if nargin < 6, q_in = []; end
    if nargin < 7, skip_above_sats = Inf; end
    result.run_idx = run_idx;
    result.worker_id = worker_id;
    result.is_candidate = false;
    result.t_faster = 0;
    result.t_fast = 0;
    result.t_detailed = 0;
    result.detailed_cov = NaN;

    Cfg.StartTime = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC') + hours((rand - 0.5) * 48);
    Cfg.SampleTime = master_config.SampleTime;
    Cfg.Min_elevation_UE = master_config.Min_elevation_UE;
    Cfg.WalkerStar = isfield(master_config, 'WalkerStar') && master_config.WalkerStar;
    Cfg.Orbit_height = local_config.Orbit_height_m;
    Cfg.Num_planes = local_config.Num_planes;
    Cfg.Sats_per_plane = local_config.Sats_per_plane;
    Cfg.Inclination = local_config.Inclination;
    Cfg.Phasing = local_config.Phasing;
    Cfg.Total_sats = local_config.Total_sats;
    Cfg.Lat_range_deg = master_config.Lat_range_deg;

    p = Cfg.Num_planes;
    s = Cfg.Sats_per_plane;
    inc = Cfg.Inclination;
    phase = Cfg.Phasing;
    t_sats = Cfg.Total_sats;
    arch_str = sprintf('%d Sats %dx%d (Inc: %.1f, Ph: %d)', t_sats, p, s, inc, phase);

    % HEARTBEAT: ULTRA
    send(q, struct('worker_id', worker_id, 'run_idx', run_idx, 'state', "ULTRA", ...
                   'is_heartbeat', true, 'arch_msg', arch_str));

    % 1. ULTRA-FAST
    t1 = tic;
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Ultrafast.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = ...
        generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Ultrafast.Num_UEs);
    m1 = constellation_simulator(Cfg, false, false, [], false);
    result.t_faster = toc(t1);
    result.faster_cov = m1.worst_coverage_percent;

    if result.faster_cov < 99
        result.t_total = result.t_faster;
        result.msg = "Failed Ultra";
        result.updated_skip_above_sats = skip_above_sats;
        return;
    end

    % HEARTBEAT: FAST
    send(q, struct('worker_id', worker_id, 'run_idx', run_idx, 'state', "FAST", 'is_heartbeat', true));
    t2 = tic;
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Fast.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = ...
        generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Fast.Num_UEs);
    m2 = constellation_simulator(Cfg, false, false, [], true); % reuse satelliteScenario handle, use matlab two body
    result.t_fast = toc(t2);

    if m2.worst_coverage_percent < 99.9
        result.t_total = result.t_faster + result.t_fast;
        result.msg = "Failed Fast";
        result.updated_skip_above_sats = skip_above_sats;
        return;
    end

    % Pre-Stage-3 check: drain q_in one last time and skip if threshold was updated
    % while Stage 2 was running (avoids starting a 200s run that can't win).
    if ~isempty(q_in)
        while true
            [msg_in, got_in] = poll(q_in, 0);
            if ~got_in, break; end
            if isfield(msg_in, 'skip_above_sats')
                skip_above_sats = min(skip_above_sats, msg_in.skip_above_sats);
            end
        end
    end
    if local_config.Total_sats > skip_above_sats
        result.t_total = result.t_faster + result.t_fast;
        result.msg = sprintf('Skipped Detailed (%d sats > threshold %d)', ...
            local_config.Total_sats, skip_above_sats);
        result.updated_skip_above_sats = skip_above_sats;
        return;
    end

    % HEARTBEAT: DETAILED
    send(q, struct('worker_id', worker_id, 'run_idx', run_idx, 'state', "DETAILED", ...
                   'is_heartbeat', true, 'arch_msg', arch_str));
    t3 = tic;
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Detailed.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = ...
        generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Detailed.Num_UEs);

    m3 = constellation_simulator(Cfg, false, false, q_in, true); % toolbox propagator (accurate, process-pool compatible)
    result.t_detailed = toc(t3);
    result.t_total = result.t_faster + result.t_fast + result.t_detailed;

    % Propagate any threshold consumed inside coverage_simulator back to our local copy
    if isfield(m3, 'consumed_threshold') && m3.consumed_threshold < skip_above_sats
        skip_above_sats = m3.consumed_threshold;
    end

    if isfield(m3, 'cancelled') && m3.cancelled
        result.msg = sprintf('Cancelled Detailed (%d sats > threshold %d)', ...
            local_config.Total_sats, skip_above_sats);
        result.updated_skip_above_sats = skip_above_sats;
        return;
    end

    result.detailed_cov = m3.worst_coverage_percent;
    result.updated_skip_above_sats = skip_above_sats;

    if m3.worst_coverage_percent > 99.999
        result.is_candidate = true;
        result.msg = sprintf("PASSED! %.4f%% | %s", m3.worst_coverage_percent, arch_str);
    else
        result.msg = sprintf("Failed Detailed %.4f%% | %s", m3.worst_coverage_percent, arch_str);
    end
end

function print_dashboard(states, runs, completed, total, wall_time, event_log, cand_count, cand_target)
    clc;
    num_w = length(states);
    cols = 6;
    rows = ceil(num_w / cols);
    fprintf('======================================================================\n');
    fprintf('  GRID SEARCH: %d / %d (%.1f%%) | Time: %.1fs | Cand: %d/%d\n', ...
        completed, total, (completed/total)*100, wall_time, cand_count, cand_target);
    fprintf('======================================================================\n');
    for r = 0:(rows-1)
        for c = 1:cols
            wid = r*cols + c;
            if wid <= num_w
                s = char(states(wid));
                s_char = s(1);
                fprintf('W%02d:[%s|%5d]  ', wid, s_char, runs(wid));
            end
        end
        fprintf('\n');
    end
    fprintf('----------------------------------------------------------------------\n');
    fprintf(' RECENT EVENTS:\n');
    for i = 1:length(event_log)
        if ~isempty(event_log{i}), fprintf(' %s\n', event_log{i}); end
    end
    fprintf('======================================================================\n');
end

function generate_profiling_report(wall_time, runs_done, n_workers, t1, t2, t3, w_counts, w_math)
    fprintf('\n================== DEEP PROFILE REPORT ==================\n');
    eff = (sum(w_math) / max(wall_time * n_workers, eps)) * 100;
    fprintf('Wall Time: %.2fs | Efficiency: %.1f%% | Throughput: %.2f r/s\n', ...
        wall_time, eff, runs_done / max(wall_time, eps));
    fprintf('Stage 1 (Ultra): Avg %.4fs | Count: %d\n', mean(t1(~isnan(t1))), sum(~isnan(t1)));
    if any(~isnan(t2)), fprintf('Stage 2 (Fast):  Avg %.4fs | Count: %d\n', mean(t2(~isnan(t2))), sum(~isnan(t2))); end
    if any(~isnan(t3)), fprintf('Stage 3 (Det):   Avg %.4fs | Count: %d\n', mean(t3(~isnan(t3))), sum(~isnan(t3))); end
    fprintf('Worker Balance: Min %d, Max %d\n', min(w_counts), max(w_counts));
    fprintf('=========================================================\n\n');
end

function write_profiling_report_to_file(report_path, base_config, wall_time, runs_done, n_workers, ...
    t1, t2, t3, w_counts, w_math, best_params, all_candidates, search_grid, status_flags, n_full_grid, n_filtered_grid)
    fid = fopen(report_path, 'w');
    if fid == -1
        warning('Could not open deep profile report file: %s', report_path);
        return;
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '================== DEEP PROFILE REPORT ==================\n');
    fprintf(fid, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    eff = (sum(w_math) / max(wall_time * n_workers, eps)) * 100;
    fprintf(fid, 'Wall Time: %.2fs\n', wall_time);
    fprintf(fid, 'Worker Efficiency: %.1f%%\n', eff);
    fprintf(fid, 'Throughput: %.2f runs/s\n', runs_done / max(wall_time, eps));

    t1_valid = t1(~isnan(t1));
    if ~isempty(t1_valid)
        fprintf(fid, 'Stage 1 (Ultra): Avg %.4fs | Count: %d\n', mean(t1_valid), numel(t1_valid));
    end
    t2_valid = t2(~isnan(t2));
    if ~isempty(t2_valid)
        fprintf(fid, 'Stage 2 (Fast):  Avg %.4fs | Count: %d\n', mean(t2_valid), numel(t2_valid));
    end
    t3_valid = t3(~isnan(t3));
    if ~isempty(t3_valid)
        fprintf(fid, 'Stage 3 (Det):   Avg %.4fs | Count: %d\n', mean(t3_valid), numel(t3_valid));
    end

    fprintf(fid, 'Worker Balance: Min %d, Max %d\n', min(w_counts), max(w_counts));
    n_s1 = numel(t1(~isnan(t1)));
    fprintf(fid, '----- Search Space Coverage -----\n');
    fprintf(fid, 'Full grid:        %d constellations\n', n_full_grid);
    fprintf(fid, 'Filtered grid:    %d constellations (>= min_sats)\n', n_filtered_grid);
    fprintf(fid, 'Stage 1 eval:     %d  (%.1f%% of filtered grid)\n', n_s1, 100*n_s1/max(n_filtered_grid,1));
    fprintf(fid, 'Early stopping:   skipped %d evals (%.1f%%)\n', ...
        n_filtered_grid - n_s1, 100*(n_filtered_grid - n_s1)/max(n_filtered_grid,1));
    fprintf(fid, 'Worker failures:  %d runs skipped\n', runs_done - n_s1); % positive if a worker died
    fprintf(fid, '=========================================================\n\n');

    fprintf(fid, '---------------- BASE CONFIG PARAMETERS ----------------\n');
    fprintf(fid, '%s\n', evalc('disp(base_config)'));

    fprintf(fid, '---------------- BEST PARAMETERS ----------------\n');
    fprintf(fid, '%s\n', regexprep(evalc('disp(best_params)'), '<.*?>', ''));

    fprintf(fid, '---------------- CANDIDATE PARAMETERS ----------------\n');
    if isempty(all_candidates)
        fprintf(fid, 'No candidates found.\n');
    else
        fprintf(fid, 'Total candidates: %d\n', height(all_candidates));
        fprintf(fid, '%s\n', regexprep(evalc('disp(all_candidates)'), '<.*?>', ''));
    end

    fprintf(fid, '---------------- EVALUATION SUMMARY ----------------\n');
    n_eval = sum(~isnan(status_flags));
    n_cand = sum(status_flags == 1 | status_flags == 2);
    fprintf(fid, 'Evaluated rows tracked: %d\n', n_eval);
    fprintf(fid, 'Rows marked candidate/minimum: %d\n', n_cand);
    fprintf(fid, 'Search grid rows: %d\n', height(search_grid));

    fprintf('Deep profile report written to: %s\n', report_path);
end

function [jx, jy] = apply_density_jitter(x, y, jitter_x_flag, jitter_y_flag)
    % Initialize with exact values
    jx = x; 
    jy = y;
    
    % Find completely identical overlapping coordinate pairs
    [~, ~, ic] = unique([x, y], 'rows');
    counts = accumarray(ic, 1);
    
    % Apply dynamic jitter based on localized density
    for i = 1:max(ic)
        n = counts(i);
        if n > 1 % Only jitter if multiple points share the EXACT same (x,y)
            idx = find(ic == i);
            
            % Dynamic spread: 1 point = 0 spread. Scales with sqrt(n) up to a max of 0.20
            spread = min(0.20, 0.02 * sqrt(n - 1)); 
            
            if jitter_x_flag
                jx(idx) = x(idx) + (rand(length(idx), 1) - 0.5) * 2 * spread;
            end
            if jitter_y_flag
                jy(idx) = y(idx) + (rand(length(idx), 1) - 0.5) * 2 * spread;
            end
        end
    end
end