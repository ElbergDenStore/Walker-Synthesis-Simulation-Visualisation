function [best_params, all_candidates] = gridsearch(master_config, orbit_height_km, plot_results, min_sats)
    %% 1. Build the Ascending Grid
    grid_data = [];
    for p = master_config.Num_Planes
        for s = master_config.Sats_Plane
            for f = 1:(p-1) % Integer Walker Phase Factor. Zero is always bad so not checked
                for inc = master_config.Inc_vec
                    grid_data = [grid_data; p, s, inc, f, p*s];
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
        stall_timeout_s = 300;
    end

    evaluated_coverage = NaN(total_runs, 1);
    detailed_coverage = NaN(total_runs, 1);
    status_flags = zeros(total_runs, 1); % 0 = Invalid, 1 = Candidate, 2 = Best
    is_completed = false(total_runs, 1);

    min_sats_found = Inf;
    candidates_found = 0;
    runs_completed = 0;

    q = parallel.pool.PollableDataQueue;
    futures = parallel.FevalFuture.empty(1, 0);

    %% 4. Submit Interleaved Chunks (Glass Cockpit model)
    fprintf('\nSubmitting interleaved chunks to %d workers...\n\n', num_workers);

    for w = 1:num_workers
        % Interleave indices so all workers start from cheaper architectures.
        indices = w:num_workers:total_runs;
        if ~isempty(indices)
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
            w_runs(wid) = data.run_idx;
            w_last_msg_s(wid) = toc(t_run_start);

            if isfield(data, 'is_worker_error') && data.is_worker_error
                error('\n[!] WORKER %d ERROR ON RUN %d: %s', wid, data.run_idx, data.error_message);
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

                % --- DETERMINISTIC STOPPING LOGIC ---
                if candidates_found >= target_num_candidates
                    sorted_candidates = sortrows(all_candidates, 'Total_sats', 'ascend');
                    worst_accepted_sats = sorted_candidates.Total_sats(target_num_candidates);

                    first_pending_idx = find(~is_completed, 1, 'first');
                    if isempty(first_pending_idx)
                        break;
                    end

                    cheapest_pending_sats = search_grid.Total_sats(first_pending_idx);
                    if cheapest_pending_sats >= worst_accepted_sats
                        fprintf('\n--- Deterministic Stop! ---\n');
                        fprintf('Top %d optimal candidates found (worst accepted sats: %d).\n', target_num_candidates, worst_accepted_sats);
                        fprintf('Cheapest pending run has %d sats (cannot improve). Canceling remaining %d runs.\n', ...
                            cheapest_pending_sats, total_runs - runs_completed);
                        cancel(futures);
                        break;
                    end
                end
            end
        end

        % Refresh dashboard every 0.5 s.
        if toc(t_last_update) > 0.5
            print_dashboard(w_states, w_runs, runs_completed, total_runs, toc(t_run_start), event_log, candidates_found, target_num_candidates);
            t_last_update = tic;
        end

        if all(strcmp({futures.State}, 'finished')), break; end

        f_states = {futures.State};
        if any(strcmp(f_states, 'finished')) || any(strcmp(f_states, 'unavailable'))
            errs = {futures.Error};
            for e = 1:length(errs)
                if strcmp(f_states{e}, 'unavailable')
                    error('\n[!] WORKER %d BECAME UNAVAILABLE (likely hard crash).', e);
                end
                if ~isempty(errs{e})
                    error('\n[!] WORKER %d CRASHED: %s', e, errs{e}.message);
                end
            end
        end

        for w = 1:length(futures)
            if strcmp(f_states{w}, 'running')
                last_s = w_last_msg_s(w);
                if last_s == 0
                    last_s = 0;
                end
                if (toc(t_run_start) - last_s) > stall_timeout_s
                    error('\n[!] WORKER %d STALLED > %.0fs (state=%s, run=%d). Set master_config.Worker_stall_timeout_s to adjust timeout.', ...
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
    out_dir = fullfile('simulation_output/gridsearch_runs', folder_name);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    report_path = fullfile(out_dir, 'deep_profile_report.txt');
    write_profiling_report_to_file(report_path, master_config, run_time, runs_completed, num_workers, ...
        t_faster_all, t_fast_all, t_detailed_all, worker_run_counts, worker_total_math_time, ...
        best_params, all_candidates, search_grid, status_flags);

    %% 8. Prepare Data For Plotting
    if plot_results
        was_evaluated = ~isnan(evaluated_coverage);
        eval_grid = search_grid(was_evaluated, :);
        eval_status = status_flags(was_evaluated); % 0=Invalid, 1=Candidate, 2=Best

        isInvalid = eval_status == 0;
        isCand = eval_status == 1;
        isBest = eval_status == 2;

        history_Loss = eval_grid.Total_sats;

        phasing_deg = (eval_grid.Phasing ./ eval_grid.Num_planes) .* 360;
        history_X = table(eval_grid.Num_planes, eval_grid.Sats_per_plane, ...
            eval_grid.Inclination, phasing_deg, ...
            'VariableNames', {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees'});

        %% Plot 1: Loss vs Inclination
        f1 = figure('Visible', 'off', 'Name', 'Sats vs Inclination', 'Color', 'w'); hold on;
        scatter(history_X.Inclination(isInvalid), history_Loss(isInvalid), 35, [0.8 0.8 0.8], 'x');
        scatter(history_X.Inclination(isCand), history_Loss(isCand), 50, [0.2 0.6 0.8], 'filled', 'MarkerEdgeColor', 'k');
        scatter(history_X.Inclination(isBest), history_Loss(isBest), 100, [1 0.8 0], 'diamond', 'filled', 'MarkerEdgeColor', 'k');
        xlabel('Inclination (deg)', 'FontWeight', 'bold'); ylabel('Num Sats', 'FontWeight', 'bold');
        title("Architecture Feasibility @ " + num2str(orbit_height_km) + " km");
        legend('Invalid', 'Candidate', 'Minimum', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f1, fullfile(out_dir, 'Inclinations_NumSats.png'), 'Resolution', 300);
        close(f1);

        %% Plot 1b: Detailed-only tradeoff (coverage vs total sats)
        was_detailed = ~isnan(detailed_coverage);
        det_grid = search_grid(was_detailed, :);
        det_cov = detailed_coverage(was_detailed);

        if ~isempty(det_cov)
            det_feasible = det_cov >= 99.999;
            det_infeasible = ~det_feasible;

            f_trade = figure('Visible', 'off', 'Name', 'Detailed Tradeoff', 'Color', 'w'); hold on;
            scatter(det_grid.Total_sats(det_infeasible), det_cov(det_infeasible), 26, [0.90 0.30 0.30], 'x', 'LineWidth', 1.0);
            scatter(det_grid.Total_sats(det_feasible), det_cov(det_feasible), 36, [0.12 0.60 0.18], 'filled', 'MarkerEdgeColor', 'k');
            xlabel('Total Satellites', 'FontWeight', 'bold');
            ylabel('Worst Coverage %', 'FontWeight', 'bold');
            title("Coverage percentage @ " + num2str(orbit_height_km) + " km");
            legend('Infeasible', 'Feasible', 'Location', 'best');
            grid on; hold off;
            exportgraphics(f_trade, fullfile(out_dir, 'Detailed_Tradeoff_Coverage_vs_Sats.png'), 'Resolution', 300);
            close(f_trade);
        end

        %% Plot 2: Architecture map (planes vs sats per plane)
        planes = history_X.Num_planes;
        sats_pp = history_X.Sats_per_plane;

        f4 = figure('Visible', 'off', 'Name', 'Architecture Map', 'Color', 'w'); hold on;
        jitter_x = planes + (rand(size(planes))-0.5)*0.4;
        jitter_y = sats_pp + (rand(size(sats_pp))-0.5)*0.4;

        scatter(jitter_x(isInvalid), jitter_y(isInvalid), 30, [0.8 0.8 0.8], 'x');
        scatter(jitter_x(isCand), jitter_y(isCand), 40, history_Loss(isCand), 'filled', 'MarkerEdgeColor', 'k');
        scatter(jitter_x(isBest), jitter_y(isBest), 90, history_Loss(isBest), 'diamond', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

        colormap('parula');
        if any(isCand) || any(isBest)
            cb = colorbar;
            cb.Label.String = 'Total Satellites';
        end
        xlabel('Num Planes', 'FontWeight', 'bold'); ylabel('Sats per Plane', 'FontWeight', 'bold');
        title("Evaluated Architectures @ " + num2str(orbit_height_km) + " km");
        legend('Invalid', 'Candidate', 'Minimum', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f4, fullfile(out_dir, 'NumPlanes_SatsPerPlane.png'), 'Resolution', 300);
        close(f4);
    end
end

% =========================================================================
% WORKER & HELPER FUNCTIONS
% =========================================================================
function evaluate_chunk(q, master_config, orbit_height_km, worker_grid, original_indices, worker_id)
    % Force workers to one thread to avoid oversubscription.
    maxNumCompThreads(1);

    for i = 1:height(worker_grid)
        local_config.Orbit_height_m = orbit_height_km * 1e3;
        local_config.Num_planes = worker_grid.Num_planes(i);
        local_config.Sats_per_plane = worker_grid.Sats_per_plane(i);
        local_config.Inclination = worker_grid.Inclination(i);
        local_config.Phasing = worker_grid.Phasing(i);
        local_config.Total_sats = worker_grid.Total_sats(i);

        try
            res = run_single_evaluation(q, local_config, master_config, original_indices(i), worker_id);
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

function result = run_single_evaluation(q, local_config, master_config, run_idx, worker_id)
    result.run_idx = run_idx;
    result.worker_id = worker_id;
    result.is_candidate = false;
    result.t_faster = 0;
    result.t_fast = 0;
    result.t_detailed = 0;
    result.detailed_cov = NaN;

    Cfg.StartTime = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
    Cfg.SampleTime = 60;
    Cfg.Min_elevation_UE = master_config.Min_elevation_UE;
    Cfg.WalkerStar = false;
    Cfg.Orbit_height = local_config.Orbit_height_m;
    Cfg.Num_planes = local_config.Num_planes;
    Cfg.Sats_per_plane = local_config.Sats_per_plane;
    Cfg.Inclination = local_config.Inclination;
    Cfg.Phasing = local_config.Phasing;
    Cfg.Total_sats = local_config.Total_sats;

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
    m1 = fast_coverage_simulator_function(Cfg, false, false, false);
    result.t_faster = toc(t1);
    result.faster_cov = m1.worst_coverage_percent;

    if result.faster_cov < 99
        result.t_total = result.t_faster;
        result.msg = "Failed Ultra";
        return;
    end

    % HEARTBEAT: FAST
    send(q, struct('worker_id', worker_id, 'run_idx', run_idx, 'state', "FAST", 'is_heartbeat', true));
    t2 = tic;
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Fast.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = ...
        generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Fast.Num_UEs);
    m2 = fast_coverage_simulator_function(Cfg, false, false, true);
    result.t_fast = toc(t2);

    if m2.worst_coverage_percent < 99.9
        result.t_total = result.t_faster + result.t_fast;
        result.msg = "Failed Fast";
        return;
    end

    % HEARTBEAT: DETAILED
    send(q, struct('worker_id', worker_id, 'run_idx', run_idx, 'state', "DETAILED", ...
                   'is_heartbeat', true, 'arch_msg', arch_str));
    t3 = tic;
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Detailed.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = ...
        generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Detailed.Num_UEs);

    m3 = coverage_simulator_function(Cfg, false, false);
    result.detailed_cov = m3.worst_coverage_percent;
    result.t_detailed = toc(t3);
    result.t_total = result.t_faster + result.t_fast + result.t_detailed;

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
    fprintf('Phase 1 (Ultra): Avg %.4fs | Count: %d\n', mean(t1(~isnan(t1))), sum(~isnan(t1)));
    if any(~isnan(t2)), fprintf('Phase 2 (Fast):  Avg %.4fs | Count: %d\n', mean(t2(~isnan(t2))), sum(~isnan(t2))); end
    if any(~isnan(t3)), fprintf('Phase 3 (Det):   Avg %.4fs | Count: %d\n', mean(t3(~isnan(t3))), sum(~isnan(t3))); end
    fprintf('Worker Balance: Min %d, Max %d\n', min(w_counts), max(w_counts));
    fprintf('=========================================================\n\n');
end

function write_profiling_report_to_file(report_path, base_config, wall_time, runs_done, n_workers, ...
    t1, t2, t3, w_counts, w_math, best_params, all_candidates, search_grid, status_flags)
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
        fprintf(fid, 'Phase 1 (Ultra): Avg %.4fs | Count: %d\n', mean(t1_valid), numel(t1_valid));
    end
    t2_valid = t2(~isnan(t2));
    if ~isempty(t2_valid)
        fprintf(fid, 'Phase 2 (Fast):  Avg %.4fs | Count: %d\n', mean(t2_valid), numel(t2_valid));
    end
    t3_valid = t3(~isnan(t3));
    if ~isempty(t3_valid)
        fprintf(fid, 'Phase 3 (Det):   Avg %.4fs | Count: %d\n', mean(t3_valid), numel(t3_valid));
    end

    fprintf(fid, 'Worker Balance: Min %d, Max %d\n', min(w_counts), max(w_counts));
    fprintf(fid, '=========================================================\n\n');

    fprintf(fid, '---------------- BASE CONFIG PARAMETERS ----------------\n');
    fprintf(fid, '%s\n', evalc('disp(base_config)'));

    fprintf(fid, '---------------- BEST PARAMETERS ----------------\n');
    fprintf(fid, '%s\n', evalc('disp(best_params)'));

    fprintf(fid, '---------------- CANDIDATE PARAMETERS ----------------\n');
    if isempty(all_candidates)
        fprintf(fid, 'No candidates found.\n');
    else
        fprintf(fid, 'Total candidates: %d\n', height(all_candidates));
        fprintf(fid, '%s\n', evalc('disp(all_candidates)'));
    end

    fprintf(fid, '---------------- EVALUATION SUMMARY ----------------\n');
    n_eval = sum(~isnan(status_flags));
    n_cand = sum(status_flags == 1 | status_flags == 2);
    fprintf(fid, 'Evaluated rows tracked: %d\n', n_eval);
    fprintf(fid, 'Rows marked candidate/minimum: %d\n', n_cand);
    fprintf(fid, 'Search grid rows: %d\n', height(search_grid));

    fprintf('Deep profile report written to: %s\n', report_path);
end