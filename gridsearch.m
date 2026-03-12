function best_params = gridsearch(Cfg,plot_results, min_sats, max_sats)
    %% 1. Build the Ascending Grid
    P_vec = 4:15; % Num Planes
    S_vec = 4:15; % Sats per Plane
    Inc_vec = linspace(70, 80, 11);
    
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
        {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Factor', 'Total_Sats'});
    
    %% 2. The Smart Filters
    % Filter out anything that isn't between min_sats and max_sats satellites
    isValidTarget = search_grid.Total_Sats >= min_sats & search_grid.Total_Sats <= max_sats;
    search_grid = search_grid(isValidTarget, :);
    
    % Sort from cheapest to most expensive!
    search_grid = sortrows(search_grid, 'Total_Sats', 'ascend');
    
    fprintf('\n=== Starting Ascending Grid Search ===\n');
    fprintf('Testing %d valid architectures between %d and %d satellites...\n\n', height(search_grid), min_sats,max_sats);
    
    %% 3. Simulate until we find the Global Minimum (Batched Parallel)
    best_params = [];
    
    % Get the current parallel pool, or create one if it doesn't exist
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool(); 
    end
    num_workers = pool.NumWorkers;
    total_runs = height(search_grid);
    
    fprintf('Using %d parallel workers for batch processing...\n', num_workers);
    found_global_minimum = false;

    % Add this right before your "for batch_start = ..." loop
    evaluated_coverage = NaN(total_runs, 1);
    % Outer loop steps forward by the number of workers
    for batch_start = 1 : num_workers : total_runs
        % Calculate where this batch ends
        batch_end = min(batch_start + num_workers - 1, total_runs);
        batch_size = batch_end - batch_start + 1;
        
        fprintf('\n--- Simulating Batch: Runs %d to %d ---\n', batch_start, batch_end);
        
        % Preallocate array to hold coverage results for this batch
        batch_coverage = zeros(batch_size, 1);
        
        % Extract just the rows for this batch
        batch_grid = search_grid(batch_start:batch_end, :);
        
        % --- THE PARALLEL LOOP ---
        parfor i = 1:batch_size
            % Local Cfg for parallel safety
            local_Cfg = Cfg;
            local_Cfg.Num_planes     = batch_grid.Num_planes(i);
            local_Cfg.Sats_per_plane = batch_grid.Sats_per_plane(i);
            local_Cfg.Inclination    = batch_grid.Inclination(i);
            local_Cfg.Phasing        = batch_grid.Phasing_Factor(i); 
            local_Cfg.Total_sats     = batch_grid.Total_Sats(i);
            
            % Run simulator
            metrics = coverage_simulator_function(local_Cfg, false, false, false);
            batch_coverage(i) = metrics.worst_coverage_percent;
            
            fprintf('Finished %dx%d (Inc: %.1f, Phase: %d) -> Cov: %.2f%%\n', ...
                local_Cfg.Num_planes, local_Cfg.Sats_per_plane, local_Cfg.Inclination, local_Cfg.Phasing, batch_coverage(i));
        end
        
        % Save batch results into the master history array
        evaluated_coverage(batch_start:batch_end) = batch_coverage;

        % Evaluate the batch results
        % valid_indices = 1; %quickly test plotting
        valid_indices = find(batch_coverage >= 99.9);
        for v = 1:length(valid_indices)
            best_local_idx = valid_indices(v);
            temp_params = batch_grid(best_local_idx, :);
            
            % --- CREATE A TEMPORARY CONFIG FOR DETAILED RUN ---
            detailed_Cfg = Cfg; 
            detailed_Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
            detailed_Cfg.StopTime   = datetime('3-Jun-2025 11:59:59', 'TimeZone', 'UTC');
            detailed_Cfg.SampleTime = 20; % seconds
            detailed_Cfg.Lat_vec = linspace(55, 85, 10); 
            detailed_Cfg.Lon_vec = linspace(-60, 30, 3);
            
            % Inject final parameters
            detailed_Cfg.Num_planes     = temp_params.Num_planes;
            detailed_Cfg.Sats_per_plane = temp_params.Sats_per_plane;
            detailed_Cfg.Inclination    = temp_params.Inclination;
            detailed_Cfg.Phasing        = temp_params.Phasing_Factor;
            detailed_Cfg.Total_sats     = temp_params.Total_Sats;
            
            fprintf('  -> High-Fidelity Test for %dx%d (Total: %d)...\n', ...
                detailed_Cfg.Num_planes, detailed_Cfg.Sats_per_plane, detailed_Cfg.Total_sats);
            
            
            detailed_metrics = coverage_simulator_function(detailed_Cfg, false, true, false); %plot_results = false; use_parallel = true; calc_link = false;
            
            if detailed_metrics.worst_coverage_percent > 99.999
                fprintf('\n======================================\n');
                fprintf('====== GLOBAL MINIMUM FOUND ==========\n');
                fprintf('======================================\n');
                fprintf('Total Satellites: %d\n', temp_params.Total_Sats);
                disp(temp_params);
            
                best_params = temp_params;
                found_global_minimum = true;
                break; 
            else
                fprintf('  -> [x] Failed high-fidelity test (Cov: %.4f%%). Moving to next candidate.\n', ...
                    detailed_metrics.worst_coverage_percent);
            end
        end
        
        % If the inner loop found the winner, break the outer batch loop too!
        if found_global_minimum
            break;
        end
    end
    if isempty(best_params)
        fprintf('\n[!] GRID SEARCH EXHAUSTED [!]\n');
        fprintf('No constellation achieved 99.9%% coverage within the satellite limit.\n');
        best_params.Num_planes     = NaN;
        best_params.Sats_per_plane = NaN;
        best_params.Inclination    = NaN;
        best_params.Phasing_Factor = NaN;
        best_params.Total_Sats     = NaN;
        return;
    end
    %% --- PREPARE DATA FOR PLOTTING ---
    if plot_results
        % 1. Extract ONLY the rows we actually simulated before breaking
        was_evaluated = ~isnan(evaluated_coverage);
        eval_grid = search_grid(was_evaluated, :);
        eval_cov = evaluated_coverage(was_evaluated);
        
        % 2. Reconstruct your specific plotting variables
        history_Loss = eval_grid.Total_Sats;
        history_Constraints = 99.9 - eval_cov;
        isValid = history_Constraints <= 0;
        
        % Rebuild history_X, converting Phasing_Factor back to Degrees
        phasing_deg = (eval_grid.Phasing_Factor ./ eval_grid.Num_planes) .* 360;
        history_X = table(eval_grid.Num_planes, eval_grid.Sats_per_plane, ...
            eval_grid.Inclination, phasing_deg, ...
            'VariableNames', {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees'});
        
        %% Create Output Directory
        date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name = sprintf('%.0f_%d_%s', Cfg.Orbit_height/1e3, best_params.Total_Sats, date_str);
        out_dir = fullfile('simulation_output/gridsearch_runs', folder_name);
        
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        
        %% --- POST-RUN VISUALIZATIONS ---
        
        % --- PLOT 1: Loss vs. Inclination ---
        f1 = figure('Visible','off','Name', 'Sats vs Inclination', 'Color', 'w'); hold on;
        scatter(history_X.Inclination(~isValid), history_Loss(~isValid), 35, [0.6 0.6 0.6], 'x', 'LineWidth', 1);
        scatter(history_X.Inclination(isValid), history_Loss(isValid), 60, history_Loss(isValid), 'filled', 'MarkerEdgeColor', 'k');
        colormap('parula'); cb = colorbar; cb.Label.String = 'Total Satellites';
        xlabel('Inclination (deg)', 'FontWeight', 'bold');
        ylabel('Num Sats', 'FontWeight', 'bold');
        title("Num Sats and Inclination @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid', 'Valid');
        grid on; hold off;
        exportgraphics(f1, fullfile(out_dir, 'Inclinations_NumSats.png'), 'Resolution', 300);
        close(f1);
        
        % --- PLOT 2: The Trade-off (Coverage vs Total Satellites) ---
        cov_history = eval_cov; 
        sat_history = history_Loss;
        
        f2 = figure('Visible','off','Name', 'Trade-off Analysis', 'Color', 'w'); hold on;
        scatter(sat_history(~isValid), cov_history(~isValid), 40, [0.8 0.3 0.3], 'x', 'LineWidth', 1.2);
        scatter(sat_history(isValid), cov_history(isValid), 60, [0.2 0.7 0.2], 'filled', 'MarkerEdgeColor', 'k');
        xlabel('Num Sats', 'FontWeight', 'bold');
        ylabel('Worst Coverage (%)', 'FontWeight', 'bold');
        title("Num Sats and Coverage @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Infeasible', 'Feasible', 'Location', 'southeast');
        grid on; hold off;
        exportgraphics(f2, fullfile(out_dir, 'NumSats_Coverage.png'), 'Resolution', 300);
        close(f2);
        
        % Extract arrays for easier plotting
        planes = history_X.Num_planes;
        sats_pp = history_X.Sats_per_plane;
        phase = history_X.Phasing_Degrees;
        inc = history_X.Inclination;
        
        % --- PLOT 4: Architecture Map (Planes vs Sats per Plane) ---
        f4 = figure('Visible','off','Name', 'Architecture Map', 'Color', 'w'); hold on;
        jitter_x = planes + (rand(size(planes))-0.5)*0.4;
        jitter_y = sats_pp + (rand(size(sats_pp))-0.5)*0.4;
        
        scatter(jitter_x(~isValid), jitter_y(~isValid), 30, [0.8 0.3 0.3], 'x');
        scatter(jitter_x(isValid), jitter_y(isValid), 70, history_Loss(isValid), 'filled', 'MarkerEdgeColor', 'k');
        
        colormap('parula'); cb = colorbar; cb.Label.String = 'Num Sats';
        xlabel('Num Planes', 'FontWeight', 'bold');
        ylabel('Sats per Plane', 'FontWeight', 'bold');
        title("Evaluated Architectures @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid', 'Valid', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f4, fullfile(out_dir, 'NumPlanes_SatsPerPlane.png'), 'Resolution', 300);
        close(f4);
        
        % --- PLOT 5: Orbital Mechanics (Phasing vs Inclination) ---
        f5 = figure('Visible','off','Name', 'Phasing vs Inclination', 'Color', 'w'); hold on;
        scatter(phase(~isValid), inc(~isValid), 30, [0.8 0.3 0.3], 'x');
        scatter(phase(isValid), inc(isValid), 70, history_Loss(isValid), 'filled', 'MarkerEdgeColor', 'k');
        colormap('parula'); cb = colorbar; cb.Label.String = 'Total Satellites';
        xlabel('Phasing (deg)', 'FontWeight', 'bold');
        ylabel('Inclination (deg)', 'FontWeight', 'bold');
        title("Phasing and Inclination @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid', 'Valid', 'Location', 'northeast');
        grid on; hold off;
        exportgraphics(f5, fullfile(out_dir, 'Phasing_vs_Inclination.png'), 'Resolution', 300);
        close(f5);
        
        % --- PLOT 6: Parallel Coordinates ---
        f6 = figure('Visible','off','Name', 'Parallel Coordinates', 'Color', 'w');
        valid_data = history_X(isValid, :);
        valid_loss = history_Loss(isValid);
        
        if height(valid_data) > 0
            valid_data.Total_Sats = valid_loss;
            coord_vars = {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees', 'Total_Sats'};
            p = parallelplot(valid_data, 'CoordinateVariables', coord_vars);
            p.Color = lines(height(valid_data));
            p.LineWidth = 4; 
            p.LineAlpha = 0.8; 
            title('Optimal Solution(s) Found');
        else
            text(0.5, 0.5, 'No valid runs to plot.', 'HorizontalAlignment', 'center', 'FontSize', 14);
            axis off;
        end
        exportgraphics(f6, fullfile(out_dir, 'Top_Solutions.png'), 'Resolution', 300);
        close(f6);
    
        %% Show high res result of best constellation
        fprintf('\nRunning detailed simulations for the optimal result...\n');
        
        % Downlink Link Budget Config FR2
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
        Cfg.SampleTime = 20; % seconds
        Cfg.Lat_vec = linspace(55, 85, 10); 
        Cfg.Lon_vec = linspace(-60, 30, 3);
        
        % Inject final parameters
        Cfg.Num_planes     = best_params.Num_planes;
        Cfg.Sats_per_plane = best_params.Sats_per_plane;
        Cfg.Inclination    = best_params.Inclination;
        Cfg.Phasing        = best_params.Phasing_Factor;
        Cfg.Total_sats     = best_params.Total_Sats;
        
        % Run the detailed simulator!
        detailed_metrics = coverage_simulator_function(Cfg, true, true, true); %plot_results = true; use_parallel = true; calc_link = true;
        
        show_interactive = false;
        save_fig = true;
        show_constellation(Cfg, show_interactive, save_fig, out_dir)
    end
end
