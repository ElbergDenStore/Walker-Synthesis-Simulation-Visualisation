function [best_params, all_candidates] = gridsearch(master_config, orbit_height_km, plot_results, min_sats)
    %% Build the Ascending Grid
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
    event_log = repmat({''}, 1, 10); % Persistent history for 10 events
    
    evaluated_coverage = NaN(total_runs, 1);
    status_flags = zeros(total_runs, 1); % 0 = Invalid, 1 = Candidate, 2 = Best
    
    min_sats_found = Inf;
    % --- SUBMIT ALL RUNS TO THE BACKGROUND ---
    fprintf('\nSubmitting %d runs to the background cluster...\n', total_runs);
    
    % Preallocate futures exactly as the documentation recommends
    futures(1:total_runs) = parallel.FevalFuture;
    
    for i = 1:total_runs
        local_config.Orbit_height_m = orbit_height_km*1e3;
        local_config.Num_planes     = search_grid.Num_planes(i);
        local_config.Sats_per_plane = search_grid.Sats_per_plane(i);
        local_config.Inclination    = search_grid.Inclination(i);
        local_config.Phasing        = search_grid.Phasing(i); 
        local_config.Total_sats     = search_grid.Total_sats(i);
        
        futures(i) = parfeval(pool, @evaluate_architecture, 1, local_config, master_config, i);
        

    end

    candidates_found = 0;
    runs_completed = 0;
    
    for i = 1:total_runs
        fprintf(sprintf("run %d",i))
        [completedIdx, result] = fetchNext(futures);

        fprintf('%s\n', result.msg);
        
        evaluated_coverage(completedIdx) = result.faster_cov;
        
        if result.is_candidate
            candidates_found = candidates_found + 1;
            temp_params = search_grid(completedIdx, :);
            temp_params.Orbit_height_km = orbit_height_km;
            status_flags(completedIdx) = 1; 
            all_candidates = [all_candidates; temp_params];
            
            % Update absolute minimum
            if temp_params.Total_sats < min_sats_found
                min_sats_found = temp_params.Total_sats;
                best_params = temp_params;
            end
            
            fprintf('*** Added to Candidates %d/%d ***\n\n', candidates_found, target_num_candidates);
            
            if candidates_found >= target_num_candidates
                % 1. Sort our current pool of successful candidates by Total_sats
                sorted_candidates = sortrows(all_candidates, 'Total_sats', 'ascend');
                
                % 2. Find the worst satellite count among our Top X targets
                worst_accepted_sats = sorted_candidates.Total_sats(target_num_candidates);
                
                % 3. Find the lowest index that has NOT finished yet
                first_pending_idx = find(~is_completed, 1);
                
                if isempty(first_pending_idx)
                    break; % Everything finished naturally
                end
                
                % 4. What is the absolute lowest satellite count that could still be returned?
                cheapest_pending_sats = search_grid.Total_sats(first_pending_idx);
                
                % 5. If the cheapest possible pending run is worse than our current targets, abort.
                if cheapest_pending_sats > worst_accepted_sats
                    fprintf('\n--- Deterministic Stop! ---\n');
                    fprintf('Top %d optimal candidates found (Max Sats: %d).\n', target_num_candidates, worst_accepted_sats);
                    fprintf('The cheapest pending run has %d sats. Canceling remaining %d runs.\n', ...
                        cheapest_pending_sats, total_runs - i);
                    
                    cancel(futures);
                    break;
                end
            end
        end
    end
    
    % --- UPGRADE BEST PARAM TO STATUS 2 ---
    best_sats_count = min(search_grid.Total_sats(status_flags == 1 | status_flags == 2));
    if ~isempty(best_sats_count)
        best_indices = find(search_grid.Total_sats == best_sats_count & status_flags > 0);
        status_flags(best_indices) = 2; 
    end

    %% Create Output Directory (for report + plots)
    date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    folder_name = sprintf('%.0f_%d_%s', Cfg.Orbit_height/1e3, best_params.Total_sats, date_str);
    out_dir = fullfile('simulation_output/gridsearch_runs', folder_name);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    % Write deep profiling report and parameter summaries to file
    report_path = fullfile(out_dir, 'deep_profile_report.txt');
    write_profiling_report_to_file(report_path, Cfg, run_time, runs_completed, num_workers, ...
        t_faster_all, t_fast_all, t_detailed_all, worker_run_counts, worker_total_math_time, ...
        best_params, all_candidates, search_grid, status_flags);

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
        folder_name = sprintf('%.0f_%d_%s', orbit_height_km, best_params.Total_sats, date_str);
        out_dir = fullfile('simulation_output/gridsearch_runs', folder_name);
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end
        
        %% --- POST-RUN VISUALIZATIONS ---
        
        % --- PLOT 1: Loss vs. Inclination ---
        f1 = figure('Visible','off','Name', 'Sats vs Inclination', 'Color', 'w'); hold on;
        scatter(history_X.Inclination(isInvalid), history_Loss(isInvalid), 35, [0.8 0.8 0.8], 'x');
        scatter(history_X.Inclination(isCand), history_Loss(isCand), 50, [0.2 0.6 0.8], 'filled', 'MarkerEdgeColor', 'k');
        scatter(history_X.Inclination(isBest), history_Loss(isBest), 100, [1 0.8 0], 'diamond', 'filled', 'MarkerEdgeColor', 'k');
        xlabel('Inclination (deg)', 'FontWeight', 'bold'); ylabel('Num Sats', 'FontWeight', 'bold');
        title("Architecture Feasibility @ " + num2str(orbit_height_km) + " km");
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
        title("Evaluated Architectures @ " + num2str(orbit_height_km) + " km");
        legend('Invalid', 'Candidate', 'Minimum', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f4, fullfile(out_dir, 'NumPlanes_SatsPerPlane.png'), 'Resolution', 300);
        close(f4);
        
        % % --- PLOT 6: Parallel Coordinates (Top 5 Candidates) ---
        % f6 = figure('Visible','off','Name', 'Parallel Coordinates', 'Color', 'w');
        % valid_mask = isCand | isBest;
        % valid_data = history_X(valid_mask, :);
        % valid_loss = history_Loss(valid_mask);
        
        % if height(valid_data) > 0
        %     valid_data.Total_sats = valid_loss;
            
        %     % 1. SORT AND SLICE: Keep only the Top 5 cheapest architectures
        %     valid_data = sortrows(valid_data, 'Total_sats', 'ascend');
        %     num_to_plot = min(5, height(valid_data));
        %     plot_data = valid_data(1:num_to_plot, :);
            
        %     % 2. Assign Status Labels for the Legend
        %     group_labels = repmat({'Candidate'}, height(plot_data), 1);
        %     min_sats = min(plot_data.Total_sats);
        %     best_indices = find(plot_data.Total_sats == min_sats);
        %     for b_idx = 1:length(best_indices)
        %         group_labels{best_indices(b_idx)} = 'Minimum';
        %     end
            
        %     % Put 'Minimum' LAST so MATLAB draws it on top
        %     plot_data.Status = categorical(group_labels, {'Candidate', 'Minimum'});
            
        %     % 3. THE INTEGER TRICK: Convert discrete columns to categorical
        %     plot_data.Num_planes      = categorical(plot_data.Num_planes);
        %     plot_data.Sats_per_plane  = categorical(plot_data.Sats_per_plane);
        %     plot_data.Inclination     = categorical(plot_data.Inclination);
        %     plot_data.Phasing_Degrees = categorical(plot_data.Phasing_Degrees);
        %     plot_data.Total_sats      = categorical(plot_data.Total_sats);
            
        %     % 4. OPTIMIZED AXIS ORDER: Total_sats at the end!
        %     coord_vars = {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees', 'Total_sats'};
        %     p = parallelplot(plot_data, 'CoordinateVariables', coord_vars, 'GroupVariable', 'Status');
            
        %     % 5. FIX THE JITTER: Turn off MATLAB's automatic spreading
        %     p.Jitter = 0;
            
        %     % 6. BULLETPROOF VISUAL HIERARCHY
        %     % Dynamically check categories to ensure Grey = Candidate and Copper = Minimum
        %     cats = categories(plot_data.Status);
        %     colors = zeros(length(cats), 3);
        %     for c_idx = 1:length(cats)
        %         if strcmp(cats{c_idx}, 'Candidate')
        %             colors(c_idx, :) = [0.75 0.75 0.75]; % Light Grey
        %         elseif strcmp(cats{c_idx}, 'Minimum')
        %             colors(c_idx, :) = [0.85 0.40 0.10]; % Bold Copper
        %         end
        %     end
            
        %     p.Color = colors;
        %     p.LineWidth = 5; 
        %     p.LineAlpha = 0.9; 
            
        %     title(sprintf('Optimal Architecture Candidates (Top %d)', num_to_plot));
        % else
        %     text(0.5, 0.5, 'No valid runs to plot.', 'HorizontalAlignment', 'center', 'FontSize', 14);
        %     axis off;
        % end
        % exportgraphics(f6, fullfile(out_dir, 'Top_Solutions.png'), 'Resolution', 300);
        % close(f6);

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
function result = evaluate_architecture(local_config, master_config, run_idx)
    result.run_idx = run_idx;
    result.is_candidate = false;
    result.t_faster = 0; result.t_fast = 0; result.t_detailed = 0;

    % Base configuration
    Cfg.StartTime       = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
    Cfg.SampleTime      = 60;
    Cfg.Min_elevation_UE = master_config.Min_elevation_UE;
    Cfg.WalkerStar      = false; % We are optimizing Walker Deltas
    Cfg.Orbit_height    = local_config.Orbit_height_m;
    Cfg.Num_planes      = local_config.Num_planes;
    Cfg.Sats_per_plane  = local_config.Sats_per_plane;
    Cfg.Inclination     = local_config.Inclination;
    Cfg.Phasing         = local_config.Phasing;
    Cfg.Total_sats      = local_config.Total_sats;
    
    p = Cfg.Num_planes; s = Cfg.Sats_per_plane;
    inc = Cfg.Inclination; phase = Cfg.Phasing;

    % Run the ultra-fast low-fidelity check
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Ultrafast.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Ultrafast.Num_UEs);
    faster_metrics = fast_coverage_simulator_function(Cfg, false, false, false);
    faster_cov = faster_metrics.worst_coverage_percent;
    result.faster_cov = faster_cov;

    if result.faster_cov < 99
        result.t_total = result.t_faster;
        result.msg = "Failed Ultra"; return; 
    end
    
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Fast.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Fast.Num_UEs);
    fast_metrics = fast_coverage_simulator_function(Cfg, false, false, true); % use SGP
    fast_cov = fast_metrics.worst_coverage_percent;
    result.fast_cov = fast_cov;

    if m2.worst_coverage_percent < 99.9
        result.t_total = result.t_faster + result.t_fast;
        result.msg = "Failed Fast"; return;
    end

    % 2. FASTer CHECKs PASSED! Run the Detailed High-Fidelity Test
    Cfg.StopTime = Cfg.StartTime + hours(master_config.Detailed.Duration_h);
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = generate_equal_ish_area_UEs(master_config.Lat_range_deg, [-180, 180], master_config.Detailed.Num_UEs);

    detailed_metrics = coverage_simulator_function(Cfg, false, false); 
    detailed_cov = detailed_metrics.worst_coverage_percent;

    % Check the final result and generate the appropriate message
    if detailed_cov > 99.999
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

function write_profiling_report_to_file(report_path, Cfg, wall_time, runs_done, n_workers, t1, t2, t3, w_counts, w_math, best_params, all_candidates, search_grid, status_flags)
    fid = fopen(report_path, 'w');
    if fid == -1
        warning('Could not open deep profile report file: %s', report_path);
        return;
    end
    cleanupObj = onCleanup(@() fclose(fid));

    fprintf(fid, '================== DEEP PROFILE REPORT ==================\n');
    fprintf(fid, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    eff = (sum(w_math) / (wall_time * n_workers)) * 100;
    fprintf(fid, 'Wall Time: %.2fs\n', wall_time);
    fprintf(fid, 'Worker Efficiency: %.1f%%\n', eff);
    fprintf(fid, 'Throughput: %.2f runs/s\n', runs_done / wall_time);

    t1_valid = t1(~isnan(t1));
    fprintf(fid, 'Phase 1 (Ultra): Avg %.4fs | Count: %d\n', mean(t1_valid), numel(t1_valid));
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
    fprintf(fid, '%s\n', evalc('disp(Cfg)'));

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