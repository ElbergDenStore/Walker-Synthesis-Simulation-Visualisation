function [best_params, all_candidates] = gridsearch(Cfg, plot_results, min_sats)
    %% 1. Build the Ascending Grid
    P_vec = 2:20; % Num Planes
    S_vec = 2:20; % Sats per Plane
    Inc_vec = linspace(70, 80, 21);
    target_num_candidates = 10;
    
    grid_data = [];
    for p = P_vec
        for s = S_vec
            for f = 1:(p-1) % Integer Walker Phase Factor
                for inc = Inc_vec
                    grid_data = [grid_data; p, s, inc, f, p*s];
                end
            end
        end
    end
    
    search_grid = array2table(grid_data, 'VariableNames', ...
        {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing', 'Total_sats'});
    
    %% 2. The Smart Filters
    isValidTarget = search_grid.Total_sats >= min_sats & search_grid.Total_sats <= min_sats * 2;
    search_grid = search_grid(isValidTarget, :);
    
    % Sort from cheapest to most expensive!
    search_grid = sortrows(search_grid, 'Total_sats', 'ascend');
    
    fprintf('\n=== Starting Ascending Grid Search ===\n');
    fprintf('Testing %d valid architectures from %d satellites...\n\n', height(search_grid), min_sats);
    
    %% 3. Setup Parallel Environment and Queue
    best_params = [];
    all_candidates = table();
    
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
    event_log = repmat({''}, 1, 10); % Persistent history for 10 events
    
    evaluated_coverage = NaN(total_runs, 1);
    status_flags = zeros(total_runs, 1); % 0 = Invalid, 1 = Candidate, 2 = Best
    
    min_sats_found = Inf;
    candidates_found = 0;
    runs_completed = 0;
    
    q = parallel.pool.PollableDataQueue;
    futures = parallel.FevalFuture.empty(1, 0);
    
    %% 4. Submit Interleaved Chunks (Memory-Sliced)
    fprintf('\nSubmitting interleaved chunks to %d workers...\n\n', num_workers);
    t_submit_start = tic;
    
    for w = 1:num_workers
        % Deal out the indices like cards so all workers evaluate cheapest first
        indices = w:num_workers:total_runs;
        if ~isempty(indices)
            % Slice grid for the worker to prevent client serialization crashes
            worker_grid = search_grid(indices, :);
            futures(w) = parfeval(pool, @evaluate_chunk, 0, q, Cfg, worker_grid, indices, w);
        end
    end
    
    %% 5. Synchronous Client Loop (The Glass Cockpit)
    t_run_start = tic;
    t_last_update = tic;

    while candidates_found < target_num_candidates
        [data, gotMsg] = poll(q, 0.01); 
        
        if gotMsg
            wid = data.worker_id;
            w_runs(wid) = data.run_idx;
            
            if isfield(data, 'is_heartbeat') && data.is_heartbeat
                % LIVE HEARTBEAT: Update status immediately
                w_states(wid) = data.state;

                % Log to permanent history if starting a 5-minute Detailed run
                if strcmp(data.state, "DETAILED")
                    new_event = sprintf('[W%02d] STARTING DETAILED: Run %d | %s', wid, data.run_idx, data.arch_msg);
                    event_log = [event_log(2:end), {new_event}];
                end 
            else
                % FINISHED RUN
                runs_completed = runs_completed + 1;
                idx = data.run_idx;
                
                % Update Profiling metrics
                t_faster_all(idx) = data.t_faster;
                if data.t_fast > 0, t_fast_all(idx) = data.t_fast; end
                if data.t_detailed > 0, t_detailed_all(idx) = data.t_detailed; end
                
                worker_run_counts(wid) = worker_run_counts(wid) + 1;
                worker_total_math_time(wid) = worker_total_math_time(wid) + data.t_total;
                
                % Reset worker to Ultra state and log coverage
                w_states(wid) = "Ultra"; 
                evaluated_coverage(idx) = data.faster_cov;
                
                % Process Successful Candidates
                if data.is_candidate
                    new_event = sprintf('[W%02d] Run %d: *** %s ***', wid, idx, data.msg);
                    event_log = [event_log(2:end), {new_event}];
                    
                    candidates_found = candidates_found + 1;
                    status_flags(idx) = 1;
                    temp_params = search_grid(idx, :);
                    all_candidates = [all_candidates; temp_params];
                    
                    if temp_params.Total_sats < min_sats_found
                        min_sats_found = temp_params.Total_sats;
                        best_params = temp_params;
                    end
                end
            end
        end

        % Refresh Dashboard every 0.5s
        if toc(t_last_update) > 0.5
            print_dashboard(w_states, w_runs, runs_completed, total_runs, toc(t_run_start), event_log, candidates_found);
            t_last_update = tic;
        end
        
        % Exit conditions
        if candidates_found >= target_num_candidates, break; end
        if all(strcmp({futures.State}, 'finished')), break; end
        
        % Fail-safe for worker crashes
        if any(strcmp({futures.State}, 'finished'))
            errs = {futures.Error};
            for e = 1:length(errs)
                if ~isempty(errs{e}), error('\n[!] WORKER %d CRASHED: %s', e, errs{e}.message); end
            end
        end
    end
    
    %% --- END OF LOOP CLEANUP ---
    cancel(futures);
    run_time = toc(t_run_start);
    
    % Print the deep profiling report before final plotting
    generate_profiling_report(run_time, runs_completed, num_workers, ...
        t_faster_all, t_fast_all, t_detailed_all, worker_run_counts, worker_total_math_time);
    
    if isempty(best_params)
        fprintf('\n[!] GRID SEARCH EXHAUSTED [!]\n');
        fprintf('No constellation achieved 99.999%% coverage within limits.\n');
        return;
    end
    
    % --- UPGRADE BEST PARAM TO STATUS 2 ---
    best_sats_count = min(search_grid.Total_sats(status_flags == 1 | status_flags == 2));
    if ~isempty(best_sats_count)
        best_indices = find(search_grid.Total_sats == best_sats_count & status_flags > 0);
        status_flags(best_indices) = 2; 
    end

    %% --- PREPARE DATA FOR PLOTTING ---
    if plot_results
        was_evaluated = ~isnan(evaluated_coverage);
        eval_grid = search_grid(was_evaluated, :);
        eval_status = status_flags(was_evaluated); % 0=Invalid, 1=Candidate, 2=Best
        
        isInvalid = eval_status == 0;
        isCand    = eval_status == 1;
        isBest    = eval_status == 2;
        
        history_Loss = eval_grid.Total_sats;
        
        phasing_deg = (eval_grid.Phasing ./ eval_grid.Num_planes) .* 360;
        history_X = table(eval_grid.Num_planes, eval_grid.Sats_per_plane, ...
            eval_grid.Inclination, phasing_deg, ...
            'VariableNames', {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees'});
        
        %% Create Output Directory
        date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name = sprintf('%.0f_%d_%s', Cfg.Orbit_height/1e3, best_params.Total_sats, date_str);
        out_dir = fullfile('simulation_output/gridsearch_runs', folder_name);
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end
        
        %% --- POST-RUN VISUALIZATIONS ---
        
        % --- PLOT 1: Loss vs. Inclination ---
        f1 = figure('Visible','off','Name', 'Sats vs Inclination', 'Color', 'w'); hold on;
        scatter(history_X.Inclination(isInvalid), history_Loss(isInvalid), 35, [0.8 0.8 0.8], 'x');
        scatter(history_X.Inclination(isCand), history_Loss(isCand), 50, [0.2 0.6 0.8], 'filled', 'MarkerEdgeColor', 'k');
        scatter(history_X.Inclination(isBest), history_Loss(isBest), 100, [1 0.8 0], 'diamond', 'filled', 'MarkerEdgeColor', 'k');
        xlabel('Inclination (deg)', 'FontWeight', 'bold'); ylabel('Num Sats', 'FontWeight', 'bold');
        title("Architecture Feasibility @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid', 'Candidate', 'Minimum', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f1, fullfile(out_dir, 'Inclinations_NumSats.png'), 'Resolution', 300);
        close(f1);
        
        % --- PLOT 4: Architecture Map (Planes vs Sats per Plane) ---
        planes = history_X.Num_planes;
        sats_pp = history_X.Sats_per_plane;
        
        f4 = figure('Visible','off','Name', 'Architecture Map', 'Color', 'w'); hold on;
        jitter_x = planes + (rand(size(planes))-0.5)*0.4;
        jitter_y = sats_pp + (rand(size(sats_pp))-0.5)*0.4;
        
        % 1. Invalid Configurations (Grey Crosses)
        scatter(jitter_x(isInvalid), jitter_y(isInvalid), 30, [0.8 0.8 0.8], 'x');
        
        % 2. Candidates (Standard Circles, color-mapped by Num Sats)
        scatter(jitter_x(isCand), jitter_y(isCand), 40, history_Loss(isCand), 'filled', 'MarkerEdgeColor', 'k');
        
        % 3. Global Minima (Large Diamonds, ALSO color-mapped, thicker border)
        scatter(jitter_x(isBest), jitter_y(isBest), 90, history_Loss(isBest), 'diamond', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        colormap('parula'); 
        if any(isCand) || any(isBest)
            cb = colorbar; cb.Label.String = 'Total Satellites';
        end
        xlabel('Num Planes', 'FontWeight', 'bold'); ylabel('Sats per Plane', 'FontWeight', 'bold');
        title("Evaluated Architectures @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid', 'Candidate', 'Minimum', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f4, fullfile(out_dir, 'NumPlanes_SatsPerPlane.png'), 'Resolution', 300);
        close(f4);
        
        % --- PLOT 6: Parallel Coordinates (Top 5 Candidates) ---
        f6 = figure('Visible','off','Name', 'Parallel Coordinates', 'Color', 'w');
        valid_mask = isCand | isBest;
        valid_data = history_X(valid_mask, :);
        valid_loss = history_Loss(valid_mask);
        
        if height(valid_data) > 0
            valid_data.Total_sats = valid_loss;
            
            % 1. SORT AND SLICE: Keep only the Top 5 cheapest architectures
            valid_data = sortrows(valid_data, 'Total_sats', 'ascend');
            num_to_plot = min(5, height(valid_data));
            plot_data = valid_data(1:num_to_plot, :);
            
            % 2. Assign Status Labels for the Legend
            group_labels = repmat({'Candidate'}, height(plot_data), 1);
            min_sats = min(plot_data.Total_sats);
            best_indices = find(plot_data.Total_sats == min_sats);
            for b_idx = 1:length(best_indices)
                group_labels{best_indices(b_idx)} = 'Minimum';
            end
            
            % Put 'Minimum' LAST so MATLAB draws it on top
            plot_data.Status = categorical(group_labels, {'Candidate', 'Minimum'});
            
            % 3. THE INTEGER TRICK: Convert discrete columns to categorical
            plot_data.Num_planes      = categorical(plot_data.Num_planes);
            plot_data.Sats_per_plane  = categorical(plot_data.Sats_per_plane);
            plot_data.Inclination     = categorical(plot_data.Inclination);
            plot_data.Phasing_Degrees = categorical(plot_data.Phasing_Degrees);
            plot_data.Total_sats      = categorical(plot_data.Total_sats);
            
            % 4. OPTIMIZED AXIS ORDER: Total_sats at the end!
            coord_vars = {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees', 'Total_sats'};
            p = parallelplot(plot_data, 'CoordinateVariables', coord_vars, 'GroupVariable', 'Status');
            
            % 5. FIX THE JITTER: Turn off MATLAB's automatic spreading
            p.Jitter = 0;
            
            % 6. BULLETPROOF VISUAL HIERARCHY
            cats = categories(plot_data.Status);
            colors = zeros(length(cats), 3);
            for c_idx = 1:length(cats)
                if strcmp(cats{c_idx}, 'Candidate')
                    colors(c_idx, :) = [0.75 0.75 0.75]; % Light Grey
                elseif strcmp(cats{c_idx}, 'Minimum')
                    colors(c_idx, :) = [0.85 0.40 0.10]; % Bold Copper
                end
            end
            
            p.Color = colors;
            p.LineWidth = 5; 
            p.LineAlpha = 0.9; 
            
            title(sprintf('Optimal Architecture Candidates (Top %d)', num_to_plot));
        else
            text(0.5, 0.5, 'No valid runs to plot.', 'HorizontalAlignment', 'center', 'FontSize', 14);
            axis off;
        end
        exportgraphics(f6, fullfile(out_dir, 'Top_Solutions.png'), 'Resolution', 300);
        close(f6);

        % %% Show high res result of best constellation
        % fprintf('\nRunning detailed Link Budget simulations for the BEST result...\n');
        % 
        % Cfg.DL.Direction = "DL";
        % Cfg.DL.Target_PFD_MHz = -115;
        % Cfg.DL.B         = 4e6;     
        % Cfg.DL.f         = 12e9;    
        % % Cfg.DL.G_tx      = 42; 
        % Cfg.DL.Tx_type   = "array";      
        % Cfg.DL.G_rx      = 30;      
        % Cfg.DL.Rx_type   = "array";
        % Cfg.DL.NF        = 5;
        % 
        % 
        % Cfg.DL.G_tx = get_adjusted_tx_gain(Cfg.Orbit_height, Cfg.Min_elevation_UE, Cfg.DL.f );
        % 
        % Cfg.DL.Max_P_tx_dBm  = PFD_calc(Cfg.DL.Target_PFD_MHz, Cfg.DL.G_tx, Cfg.DL.B, Cfg.Orbit_height, Cfg.Min_elevation_UE);
        % Cfg.DL.Max_EIRP_dBm  = Cfg.DL.Max_P_tx_dBm + Cfg.DL.G_tx;
        % 
        % Cfg.Save_dir = out_dir;
        % Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
        % Cfg.StopTime   = datetime('3-Jun-2025 11:59:59', 'TimeZone', 'UTC');
        % Cfg.SampleTime = 20;
        % Cfg.Lat_vec = linspace(55, 85, 10); 
        % Cfg.Lon_vec = linspace(-60, 30, 3);
        % Cfg.Equal_UE_area = true;
        % Cfg.FRF = 3;
        % Cfg.RU = 1;
        % 
        % Cfg.Num_planes     = best_params.Num_planes;
        % Cfg.Sats_per_plane = best_params.Sats_per_plane;
        % Cfg.Inclination    = best_params.Inclination;
        % Cfg.Phasing        = best_params.Phasing;
        % Cfg.Total_sats     = best_params.Total_sats;
        % 
        % % Run the detailed simulator!
        % detailed_metrics = coverage_simulator_function(Cfg, true, true, true); %plot_results = true; use_parallel = true; calc_link = true;
        % 
        % show_constellation(Cfg, false, true, out_dir);%Show interactive = false; savefig = true;
    end
end

% =========================================================================
% WORKER & HELPER FUNCTIONS
% =========================================================================
function evaluate_chunk(q, base_Cfg, worker_grid, original_indices, worker_id)
    % FORCE single-threading so workers don't clash on CPU resources
    maxNumCompThreads(1); 
    
    for i = 1:height(worker_grid)
        local_Cfg = base_Cfg;
        local_Cfg.Num_planes     = worker_grid.Num_planes(i);
        local_Cfg.Sats_per_plane = worker_grid.Sats_per_plane(i);
        local_Cfg.Inclination    = worker_grid.Inclination(i);
        local_Cfg.Phasing        = worker_grid.Phasing(i); 
        local_Cfg.Total_sats     = worker_grid.Total_sats(i);
        
        res = run_single_evaluation(q, local_Cfg, original_indices(i), worker_id);
        res.is_heartbeat = false;
        send(q, res);
    end
end

function result = run_single_evaluation(q, Cfg, run_idx, worker_id)
    result.run_idx = run_idx; result.worker_id = worker_id;
    result.is_candidate = false;
    result.t_faster = 0; result.t_fast = 0; result.t_detailed = 0;

    p = Cfg.Num_planes; s = Cfg.Sats_per_plane;
    inc = Cfg.Inclination; phase = Cfg.Phasing;
    t_sats = Cfg.Total_sats;
    arch_str = sprintf('%d Sats %dx%d (Inc: %.1f, Ph: %d)', t_sats, p, s, inc, phase);
    
    % 1. ULTRA-FAST
    t1 = tic;
    Cfg.StopTime = Cfg.StartTime + hours(2);
    m1 = fast_coverage_simulator_function(Cfg, false, false, false);
    result.t_faster = toc(t1);
    result.faster_cov = m1.worst_coverage_percent;

    if result.faster_cov < 99
        result.t_total = result.t_faster;
        result.msg = "Failed Ultra"; return; 
    end
    
    % HEARTBEAT: FAST
    send(q, struct('worker_id', worker_id, 'run_idx', run_idx, 'state', "FAST", 'is_heartbeat', true));
    t2 = tic;
    Cfg.StopTime = Cfg.StartTime + hours(24);
    m2 = fast_coverage_simulator_function(Cfg, false, false, true); 
    result.t_fast = toc(t2);

    if m2.worst_coverage_percent < 99.9
        result.t_total = result.t_faster + result.t_fast;
        result.msg = "Failed Fast"; return;
    end

    % HEARTBEAT: DETAILED
    send(q, struct('worker_id', worker_id, 'run_idx', run_idx, 'state', "DETAILED", ...
                   'is_heartbeat', true, 'arch_msg', arch_str));
    t3 = tic;
    % [Detailed Config Setup]
    detailed_Cfg = Cfg; 
    detailed_Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
    detailed_Cfg.StopTime   = datetime('2-Jun-2025 11:59:59', 'TimeZone', 'UTC');
    detailed_Cfg.SampleTime = 60; 
    detailed_Cfg.Lat_vec = linspace(55, 85, 14); 
    detailed_Cfg.Lon_vec = linspace(-180, 180, 3);
    detailed_Cfg.Equal_UE_area = true; 
    
    m3 = coverage_simulator_function(detailed_Cfg, false, false); 
    result.t_detailed = toc(t3);
    result.t_total = result.t_faster + result.t_fast + result.t_detailed;
    
    if m3.worst_coverage_percent > 99.999
        result.is_candidate = true;
        result.msg = sprintf("PASSED! %.4f%% | %s", m3.worst_coverage_percent, arch_str);
    else
        result.msg = sprintf("Failed Detailed %.4f%% | %s", m3.worst_coverage_percent, arch_str);
    end
end

function print_dashboard(states, runs, completed, total, wall_time, event_log, cand_count)
    clc;
    num_w = length(states); cols = 6; rows = ceil(num_w / cols);
    fprintf('======================================================================\n');
    fprintf('  GRID SEARCH: %d / %d (%.1f%%) | Time: %.1fs | Cand: %d/10\n', ...
        completed, total, (completed/total)*100, wall_time, cand_count);
    fprintf('======================================================================\n');
    for r = 0:(rows-1)
        for c = 1:cols
            wid = r*cols + c;
            if wid <= num_w
                s = char(states(wid)); s_char = s(1);
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
    eff = (sum(w_math) / (wall_time * n_workers)) * 100;
    fprintf('Wall Time: %.2fs | Efficiency: %.1f%% | Throughput: %.2f r/s\n', wall_time, eff, runs_done/wall_time);
    fprintf('Phase 1 (Ultra): Avg %.4fs | Count: %d\n', mean(t1(~isnan(t1))), sum(~isnan(t1)));
    if any(~isnan(t2)), fprintf('Phase 2 (Fast):  Avg %.4fs | Count: %d\n', mean(t2(~isnan(t2))), sum(~isnan(t2))); end
    if any(~isnan(t3)), fprintf('Phase 3 (Det):   Avg %.4fs | Count: %d\n', mean(t3(~isnan(t3))), sum(~isnan(t3))); end
    fprintf('Worker Balance: Min %d, Max %d\n', min(w_counts), max(w_counts));
    fprintf('=========================================================\n\n');
end