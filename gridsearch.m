function [best_params, all_candidates] = gridsearch(Cfg, plot_results, min_sats)
    %% 1. Force kill the current pool
    delete(gcp('nocreate')); % necessary or it will get stuck

    %% 1. Build the Ascending Grid
    P_vec = 4:15; % Num Planes % 15 both places makes sense to me
    S_vec = 4:15; % Sats per Plane
    Inc_vec = linspace(70, 80, 11);
    target_num_candidates = 10;
    
    grid_data = [];
    for p = P_vec
        for s = S_vec
            for f = 1:(p-1) % Integer Walker Phase Factor. Zero is always bad so not checked
                for inc = Inc_vec
                    grid_data = [grid_data; p, s, inc, f, p*s];
                end
            end
        end
    end
    
    search_grid = array2table(grid_data, 'VariableNames', ...
        {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing', 'Total_sats'});
    
    %% 2. The Smart Filters
    isValidTarget = search_grid.Total_sats >= min_sats;
    search_grid = search_grid(isValidTarget, :);
    
    % Sort from cheapest to most expensive!
    search_grid = sortrows(search_grid, 'Total_sats', 'ascend');
    
    fprintf('\n=== Starting Ascending Grid Search ===\n');
    fprintf('Testing %d valid architectures from %d satellites...\n\n', height(search_grid), min_sats);
    
    %% 3. Simulate until we find the Global Minimum + N candidates
    best_params = [];
    all_candidates = table();
    
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool(); end
    num_workers = pool.NumWorkers;
    total_runs = height(search_grid);
    
    fprintf('Using %d parallel workers for batch processing...\n', num_workers);
    
    evaluated_coverage = NaN(total_runs, 1);
    status_flags = zeros(total_runs, 1); % 0 = Invalid, 1 = Candidate, 2 = Best
    
    min_sats_found = Inf;
    % --- SUBMIT ALL RUNS TO THE BACKGROUND ---
    fprintf('\nSubmitting %d runs to the background cluster...\n', total_runs);
    
    % Preallocate futures exactly as the documentation recommends
    futures(1:total_runs) = parallel.FevalFuture;
    
    % Add this safety net
    % cleanupObj = onCleanup(@() cancel(futures)); %does not work anyway
    
    for i = 1:total_runs
        local_Cfg = Cfg;
        local_Cfg.Num_planes     = search_grid.Num_planes(i);
        local_Cfg.Sats_per_plane = search_grid.Sats_per_plane(i);
        local_Cfg.Inclination    = search_grid.Inclination(i);
        local_Cfg.Phasing        = search_grid.Phasing(i); 
        local_Cfg.Total_sats     = search_grid.Total_sats(i);
        
        % Notice: We deleted the 'dq' input!
        futures(i) = parfeval(pool, @evaluate_architecture, 1, local_Cfg, i);
        
        % Keep the UI alive during submission
        % if mod(i, 50) == 0
        %     pause(0.01);
        %     fprintf(sprintf("%d \n",i))
        % end
    end
    % fprintf("done with le futures")
    
    % --- FETCH RESULTS LIVE AS THEY FINISH ---
    candidates_found = 0;
    
    % fetchNext waits until ANY worker finishes, then hands us the result instantly!
    for i = 1:total_runs
        fprintf(sprintf("run %d",i))
        try
            [completedIdx, result] = fetchNext(futures);
            
            % 1. PRINT THE WORKER's MESSAGE LIVE!
            fprintf('%s\n', result.msg);
            
            % 2. Store the data
            evaluated_coverage(completedIdx) = result.faster_cov;
            
            if result.is_candidate
                candidates_found = candidates_found + 1;
                temp_params = search_grid(completedIdx, :);
                status_flags(completedIdx) = 1; 
                all_candidates = [all_candidates; temp_params];
                
                % Update absolute minimum
                if temp_params.Total_sats < min_sats_found
                    min_sats_found = temp_params.Total_sats;
                    best_params = temp_params;
                end
                
                fprintf('*** Added to Candidates %d/%d ***\n\n', candidates_found, target_num_candidates);
                
                % 3. ABORT IF WE HIT THE TARGET
                if candidates_found >= target_num_candidates
                    fprintf('\n--- Target of %d candidates reached! Canceling remaining runs. ---\n', target_num_candidates);
                    cancel(futures);
                    break;
                end
            end
            
        catch ME
            % If a worker crashes inside, fetchNext throws an error. Catch and print it live!
            fprintf('\n[!] CRASH IN BACKGROUND RUN: %s\n', ME.message);
        end
    end
    
    if isempty(best_params)
        fprintf('\n[!] GRID SEARCH EXHAUSTED [!]\n');
        fprintf('No constellation achieved 99.999%% coverage within limits.\n');
        return;
    end
    
    % --- UPGRADE BEST PARAM TO STATUS 2 ---
    % Find ALL configurations that share the absolute minimum satellite count
    best_sats_count = min(search_grid.Total_sats(status_flags == 1 | status_flags == 2));
    if ~isempty(best_sats_count)
        best_indices = find(search_grid.Total_sats == best_sats_count & status_flags > 0);
        status_flags(best_indices) = 2; % Mark ALL of them as the global minimum
    end

    %% --- PREPARE DATA FOR PLOTTING ---
    if plot_results
        was_evaluated = ~isnan(evaluated_coverage);
        eval_grid = search_grid(was_evaluated, :);
        eval_cov = evaluated_coverage(was_evaluated);
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
            % Dynamically check categories to ensure Grey = Candidate and Copper = Minimum
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
% WORKER FUNCTION (Executes invisibly, returns a single result struct)
% =========================================================================
function result = evaluate_architecture(Cfg, run_idx)
    result.run_idx = run_idx;
    result.is_candidate = false;

    p = Cfg.Num_planes; s = Cfg.Sats_per_plane;
    inc = Cfg.Inclination; phase = Cfg.Phasing;

    % 1. Run the ultra-fast low-fidelity check
    Cfg.StopTime = Cfg.StartTime + hours(2);
    faster_metrics = fast_coverage_simulator_function(Cfg, false, false, false, false);
    faster_cov = faster_metrics.worst_coverage_percent;
    result.faster_cov = faster_cov;

    % If it fails, write the failure message and return immediately
    if faster_cov < 97
        result.msg = sprintf('[-] %dx%d (Inc: %.1f, Phase: %d) -> Failed Ultra Fast (Cov: %.2f%%)', p, s, inc, phase, faster_cov);
        return;
    end
    
    Cfg.StopTime = Cfg.StartTime + hours(24);
    fast_metrics = fast_coverage_simulator_function(Cfg, false, false, false, true); % use SGP
    fast_cov = fast_metrics.worst_coverage_percent;
    result.fast_cov = fast_cov;

    % If it fails again, write the failure message and return immediately
    if fast_cov < 99.9
        result.msg = sprintf('[-] %dx%d (Inc: %.1f, Phase: %d) -> Failed Fast (Cov: %.2f%%)', p, s, inc, phase, fast_cov);
        result.msg = sprintf('[-] %dx%d (Inc: %.1f, Phase: %d) -> Passed Ultra Fast (%.2f%%), Failed Fast (%.4f%%)', p, s, inc, phase, faster_cov, fast_cov);
        return;
    end

    % 2. FASTer CHECKs PASSED! Run the Detailed High-Fidelity Test
    detailed_Cfg = Cfg; 
    detailed_Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
    detailed_Cfg.StopTime   = datetime('3-Jun-2025 11:59:59', 'TimeZone', 'UTC');
    detailed_Cfg.SampleTime = 20; 
    detailed_Cfg.Lat_vec = linspace(55, 85, 10); 
    detailed_Cfg.Lon_vec = linspace(-60, 30, 3);
    detailed_Cfg.Equal_UE_area = true; 

    detailed_metrics = coverage_simulator_function(detailed_Cfg, false, false, false); 
    detailed_cov = detailed_metrics.worst_coverage_percent;

    % Check the final result and generate the appropriate message
    if detailed_cov > 99.999
        result.is_candidate = true;
        result.msg = sprintf('[+] %dx%d (Inc: %.1f, Phase: %d) -> PASSED BOTH! (Detailed Cov: %.4f%%)', p, s, inc, phase, detailed_cov);
    else
        result.msg = sprintf('[-] %dx%d (Inc: %.1f, Phase: %d) -> Passed Fast (%.2f%%), Failed Detailed (%.4f%%)', p, s, inc, phase, fast_cov, detailed_cov);
    end
end
