function [best_params, all_candidates] = gridsearch(Cfg, plot_results, min_sats, max_sats)
    %% 1. Build the Ascending Grid
    P_vec = 4:15; % Num Planes
    S_vec = 4:15; % Sats per Plane
    Inc_vec = linspace(70, 80, 11);
    target_num_candidates = 10;
    % extra_search_percent = 15; % percent extra to look for solutions
    
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
        {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Factor', 'Total_Sats'});
    
    %% 2. The Smart Filters
    isValidTarget = search_grid.Total_Sats >= min_sats & search_grid.Total_Sats <= max_sats;
    search_grid = search_grid(isValidTarget, :);
    
    % Sort from cheapest to most expensive!
    search_grid = sortrows(search_grid, 'Total_Sats', 'ascend');
    
    fprintf('\n=== Starting Ascending Grid Search ===\n');
    fprintf('Testing %d valid architectures between %d and %d satellites...\n\n', height(search_grid), min_sats,max_sats);
    
    %% 3. Simulate until we find the Global Minimum + 5% Margin (Batched Parallel)
    best_params = [];
    all_candidates = table(); % NEW: Store all viable options
    
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool(); end
    num_workers = pool.NumWorkers;
    total_runs = height(search_grid);
    
    fprintf('Using %d parallel workers for batch processing...\n', num_workers);
    
    evaluated_coverage = NaN(total_runs, 1);
    status_flags = zeros(total_runs, 1); % 0 = Invalid, 1 = Candidate, 2 = Best
    
    min_sats_found = Inf;
    max_sats_to_check = Inf;
    
    for batch_start = 1 : num_workers : total_runs
        % --- NEW BREAK CONDITION ---
        % Because the grid is sorted by Total_Sats, if the current batch starts 
        % higher than our 5% limit, we know we're done!
        % if search_grid.Total_Sats(batch_start) > max_sats_to_check
        %     fprintf('\n--- Exceeded %d\% limit above minimum found (%d sats). Stopping search! ---\n',extra_search_percent, max_sats_to_check);
        %     break;
        % end
        if height(all_candidates) >= target_num_candidates
            fprintf('\n--- Found %d viable candidates. Target reached! Stopping search. ---\n', height(all_candidates));
            break;
        end
        
        batch_end = min(batch_start + num_workers - 1, total_runs);
        batch_size = batch_end - batch_start + 1;
        
        fprintf('\n--- Simulating Batch: Runs %d to %d ---\n', batch_start, batch_end);
        
        batch_coverage = zeros(batch_size, 1);
        batch_grid = search_grid(batch_start:batch_end, :);
        
        % --- THE PARALLEL LOOP ---
        parfor i = 1:batch_size
            local_Cfg = Cfg;
            local_Cfg.Num_planes     = batch_grid.Num_planes(i);
            local_Cfg.Sats_per_plane = batch_grid.Sats_per_plane(i);
            local_Cfg.Inclination    = batch_grid.Inclination(i);
            local_Cfg.Phasing        = batch_grid.Phasing_Factor(i); 
            local_Cfg.Total_sats     = batch_grid.Total_Sats(i);
            
            metrics = coverage_simulator_function(local_Cfg, false, false, false);
            batch_coverage(i) = metrics.worst_coverage_percent;
            
            fprintf('Finished %dx%d (Inc: %.1f, Phase: %d) -> Cov: %.2f%%\n', ...
                local_Cfg.Num_planes, local_Cfg.Sats_per_plane, local_Cfg.Inclination, local_Cfg.Phasing, batch_coverage(i));
        end
        
        evaluated_coverage(batch_start:batch_end) = batch_coverage;
        valid_indices = find(batch_coverage >= 99.9);
        
        for v = 1:length(valid_indices)
            best_local_idx = valid_indices(v);
            temp_params = batch_grid(best_local_idx, :);
            global_row_idx = batch_start + best_local_idx - 1;
            
            % --- DETAILED HIGH FIDELITY RUN ---
            detailed_Cfg = Cfg; 
            detailed_Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
            detailed_Cfg.StopTime   = datetime('3-Jun-2025 11:59:59', 'TimeZone', 'UTC');
            detailed_Cfg.SampleTime = 20; 
            detailed_Cfg.Lat_vec = linspace(55, 85, 10); 
            detailed_Cfg.Lon_vec = linspace(-60, 30, 3);
            
            detailed_Cfg.Num_planes     = temp_params.Num_planes;
            detailed_Cfg.Sats_per_plane = temp_params.Sats_per_plane;
            detailed_Cfg.Inclination    = temp_params.Inclination;
            detailed_Cfg.Phasing        = temp_params.Phasing_Factor;
            detailed_Cfg.Total_sats     = temp_params.Total_Sats;
            
            fprintf('  -> High-Fidelity Test for %dx%d (Total: %d)...\n', ...
                detailed_Cfg.Num_planes, detailed_Cfg.Sats_per_plane, detailed_Cfg.Total_sats);
            
            detailed_metrics = coverage_simulator_function(detailed_Cfg, false, true, false); %plot_results = false; use_parallel = true; calc_link = false;
            
            if detailed_metrics.worst_coverage_percent > 99.999
                fprintf('  -> PASSED. Added to Candidates.\n');
                
                % Mark as a candidate
                status_flags(global_row_idx) = 1; 
                all_candidates = [all_candidates; temp_params];
                
                % Is this the new absolute best/cheapest?
                if temp_params.Total_Sats < min_sats_found
                    min_sats_found = temp_params.Total_Sats;
                    best_params = temp_params;
                    
                    fprintf('\n======================================================\n');
                    fprintf('New Global Minimum for %d km Found: %d Sats.\n', (Cfg.Orbit_height / 1000), min_sats_found);
                    fprintf('Current Candidates Pool: %d / %d\n', height(all_candidates), target_num_candidates);
                    fprintf('======================================================\n');
                end
            else
                fprintf('  -> Failed high-fidelity test (Cov: %.4f%%).\n', detailed_metrics.worst_coverage_percent);
            end
        end
    end
    
    if isempty(best_params)
        fprintf('\n[!] GRID SEARCH EXHAUSTED [!]\n');
        fprintf('No constellation achieved 99.999%% coverage within limits.\n');
        return;
    end
    
    % --- UPGRADE BEST PARAM TO STATUS 2 ---
    % Find ALL configurations that share the absolute minimum satellite count
    best_sats_count = min(search_grid.Total_Sats(status_flags == 1 | status_flags == 2));
    if ~isempty(best_sats_count)
        best_indices = find(search_grid.Total_Sats == best_sats_count & status_flags > 0);
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
        
        history_Loss = eval_grid.Total_Sats;
        
        phasing_deg = (eval_grid.Phasing_Factor ./ eval_grid.Num_planes) .* 360;
        history_X = table(eval_grid.Num_planes, eval_grid.Sats_per_plane, ...
            eval_grid.Inclination, phasing_deg, ...
            'VariableNames', {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees'});
        
        %% Create Output Directory
        date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name = sprintf('%.0f_%d_%s', Cfg.Orbit_height/1e3, best_params.Total_Sats, date_str);
        out_dir = fullfile('simulation_output/gridsearch_runs', folder_name);
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end
        
        %% --- POST-RUN VISUALIZATIONS ---
        
        % --- PLOT 1: Loss vs. Inclination ---
        f1 = figure('Visible','off','Name', 'Sats vs Inclination', 'Color', 'w'); hold on;
        scatter(history_X.Inclination(isInvalid), history_Loss(isInvalid), 35, [0.8 0.8 0.8], 'x');
        scatter(history_X.Inclination(isCand), history_Loss(isCand), 50, [0.2 0.6 0.8], 'filled', 'MarkerEdgeColor', 'k');
        scatter(history_X.Inclination(isBest), history_Loss(isBest), 200, [1 0.8 0], 'pentagram', 'filled', 'MarkerEdgeColor', 'k');
        xlabel('Inclination (deg)', 'FontWeight', 'bold'); ylabel('Num Sats', 'FontWeight', 'bold');
        title("Architecture Feasibility @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid', 'Candidate', 'Global Minimum', 'Location', 'best');
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
        legend('Invalid', 'Candidate', 'Global Minimum', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f4, fullfile(out_dir, 'NumPlanes_SatsPerPlane.png'), 'Resolution', 300);
        close(f4);
        
        % --- PLOT 6: Parallel Coordinates (Candidates Only) ---
        f6 = figure('Visible','off','Name', 'Parallel Coordinates', 'Color', 'w');
        valid_mask = isCand | isBest;
        valid_data = history_X(valid_mask, :);
        valid_loss = history_Loss(valid_mask);
        
        if height(valid_data) > 0
            valid_data.Total_Sats = valid_loss;
            
            % 1. Create a grouping array for the legend and colors
            group_labels = repmat({'Candidate'}, height(valid_data), 1);
            
            % 2. Find the absolute minimum and overwrite their labels
            min_sats = min(valid_loss);
            best_local_indices = find(valid_loss == min_sats);
            for b_idx = 1:length(best_local_indices)
                group_labels{best_local_indices(b_idx)} = 'Global Minimum';
            end
            
            % 3. Add to the table as a categorical variable
            valid_data.Status = categorical(group_labels);
            
            % Force the order so we know exactly which color goes to which group
            valid_data.Status = reordercats(valid_data.Status, {'Candidate', 'Global Minimum'});
            
            % 4. Plot using the GroupVariable
            coord_vars = {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees', 'Total_Sats'};
            p = parallelplot(valid_data, 'CoordinateVariables', coord_vars, 'GroupVariable', 'Status');
            
            % 5. Apply the two distinct colors mapping to the two categories
            % Row 1: Muted Slate Blue (Candidate)
            % Row 2: Copper/Orange (Global Minimum)
            p.Color = [0.4 0.5 0.6; 0.85 0.40 0.10]; 
            p.LineWidth = 3; 
            p.LineAlpha = 0.8; 
            
            title('Optimal Candidates Found (<= 5% of Minimum)');
        else
            text(0.5, 0.5, 'No valid runs to plot.', 'HorizontalAlignment', 'center', 'FontSize', 14);
            axis off;
        end
        exportgraphics(f6, fullfile(out_dir, 'Top_Solutions.png'), 'Resolution', 300);
        close(f6);

        %% Show high res result of best constellation
        fprintf('\nRunning detailed Link Budget simulations for the BEST result...\n');
        
        Cfg.DL.Direction = "DL";
        Cfg.DL.B         = 2e6;     
        Cfg.DL.f         = 20e9;    
        Cfg.DL.P_tx_dBm  = 20; 
        Cfg.DL.G_tx      = 42; 
        Cfg.DL.Tx_type   = "array";      
        Cfg.DL.G_rx      = 32;      
        Cfg.DL.Rx_type   = "array";      
        Cfg.DL.NF        = 5;
        Cfg.DL.EIRP_dBm  = Cfg.DL.P_tx_dBm + Cfg.DL.G_tx;
        
        Cfg.Save_dir = out_dir;
        Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
        Cfg.StopTime   = datetime('3-Jun-2025 11:59:59', 'TimeZone', 'UTC');
        Cfg.SampleTime = 20;
        Cfg.Lat_vec = linspace(55, 85, 10); 
        Cfg.Lon_vec = linspace(-60, 30, 3);
        
        Cfg.Num_planes     = best_params.Num_planes;
        Cfg.Sats_per_plane = best_params.Sats_per_plane;
        Cfg.Inclination    = best_params.Inclination;
        Cfg.Phasing        = best_params.Phasing_Factor;
        Cfg.Total_sats     = best_params.Total_Sats;
        
        % Run the detailed simulator!
        detailed_metrics = coverage_simulator_function(Cfg, true, true, true); %plot_results = true; use_parallel = true; calc_link = true;
        
        show_constellation(Cfg, false, true, out_dir);%Show interactive = false; savefig = true;
        
    end
end