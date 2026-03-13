function best_params = surrogateopt_constellation_optimizer(orbit_height, num_runs, plot_results)

    Cfg.StartTime  = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
    Cfg.StopTime   = datetime('1-Jun-2025 11:59:59', 'TimeZone', 'UTC'); % Runtime linearly dependent on sim time and sampletime
    Cfg.SampleTime = 60; % seconds
    Cfg.Lat_vec = linspace(55, 85, 6); % Runtime dependent on lat times lon
    Cfg.Lon_vec = linspace(-60, 30, 1);
    Cfg.Min_elevation_UE = 20;
    Cfg.WalkerStar     = false;
    Cfg.Orbit_height = orbit_height;
    
    %% Define the search space
    var_planes = optimizableVariable('Num_planes', [4, 15], 'Type', 'integer');
    var_sats   = optimizableVariable('Sats_per_plane', [4, 15], 'Type', 'integer');
    var_inc    = optimizableVariable('Inclination', [70, 80]);
    var_phasing = optimizableVariable('Phasing_Degrees', [0, 359]);
    
    % search_vars = [var_planes, var_sats, var_inc, var_phasing];

    %% Define the search space (Vectorized for surrogateopt)
    % Variables: [Num_planes, Sats_per_plane, Inclination, Phasing_Factor]
    % Note: Phasing_Factor is an integer from 0 to 14 (max planes - 1)
    % Lower bounds: Min 4 planes, min 4 sats/plane, min 70 inc, min 0 phase
    lb = [4, 4, 70, 0]; 
    
    % Upper bounds: Max 15 planes, max 15 sats/plane, max 80 inc, max 14 phase
    ub = [15, 15, 80, 360]; 
    
    % Integer constraints: Variables 1, 2 CANNOT be decimals
    intcon =[1, 2];

    options = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', num_runs, ... 
        'UseParallel', true, ...
        'PlotFcn', 'surrogateoptplot');     % This UI has a real Stop button
    
    % surrogateopt passes a numeric array, so we wrap it for your simulator
    obj_fun_surrogate = @(x) surrogateWrapper(x, Cfg); 
    
    % Run optimizer
    [best_x, best_fval, exitflag, output, trials] = surrogateopt(obj_fun_surrogate, lb, ub, intcon, options);
    
    %% Extract and Save the Best Results
    % Convert the best numeric array back into a table so your output matches
    best_params = table(best_x(1), best_x(2), best_x(3), best_x(4), ...
        'VariableNames', {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Factor'});
    
    fprintf('\n=== OPTIMIZATION COMPLETE ===\n');
    fprintf('Best Constellation Found:\n');
    disp(best_params);
    
    %% --- TRANSLATE RESULTS FOR YOUR EXISTING PLOTS ---
    valid_evals = ~isnan(trials.Fval);
    history_Loss = trials.Fval(valid_evals);
    history_Constraints = trials.Ineq(valid_evals);
    
    X_mat = trials.X(valid_evals, :);
    
    % We convert the integer Phasing_Factor back into continuous degrees 
    % JUST so your existing plotting code doesn't break!
    planes_array = X_mat(:, 1);
    phasing_factor_array = mod(X_mat(:, 4), planes_array);
    phasing_degrees_array = (phasing_factor_array ./ planes_array) .* 360;
    
    % Reconstruct history_X table exactly as your plotting code expects it
    history_X = table(X_mat(:,1), X_mat(:,2), X_mat(:,3), phasing_degrees_array, ...
        'VariableNames', {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees'});
    
    cov_history = 99.9 - history_Constraints;
    sat_history = history_Loss;
    
    best_idx = find(history_Loss == best_fval & history_Constraints <= 0, 1, 'first');
    fprintf('--- Metrics for the Best Valid Constellation ---\n');
    fprintf('Total Sats: %d\nWorst Coverage: %.2f%%\n', sat_history(best_idx), cov_history(best_idx))

     if plot_results
        %% Create Output Directory
        % Format: simulation_output/bayesian_runs/orbit_sats_inclination_phasing_date
        date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name = sprintf('%.0f_%d_%d_%s', ...
            Cfg.Orbit_height/1e3, best_params.Num_planes,best_params.Sats_per_plane, date_str);
        out_dir = fullfile('simulation_output/bayesian_runs', folder_name);
        
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end

        %% --- POST-RUN VISUALIZATIONS ---
        % history_X = results.XTrace; 
        % history_Loss = results.ObjectiveTrace;
        % history_Constraints = results.ConstraintsTrace;
        
        % Mask for valid runs (99.9% coverage or better)
        isValid = history_Constraints <= 0; 
        
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
        % cov_history = cellfun(@(x) x.Cov_percent, results.UserDataTrace);
        % sat_history = cellfun(@(x) x.Total_sats, results.UserDataTrace);
        
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

        % --- PLOT 3: Best Valid Loss Over Time ---
        % Calculate the running minimum of valid objective evaluations
        num_evals = length(history_Loss);
        best_valid_so_far = NaN(num_evals, 1);
        current_best = inf;
        
        for k = 1:num_evals
            % Update current best only if the run is valid AND lower than previous best
            if history_Constraints(k) <= 0 && history_Loss(k) < current_best
                current_best = history_Loss(k);
            end
            
            % Record the best valid loss found up to run 'k'
            if ~isinf(current_best)
                best_valid_so_far(k) = current_best;
            end
        end
        
        f3 = figure('Visible','off','Name', 'Best Valid Loss Over Time', 'Color', 'w'); hold on;
        
        % Plot all runs in the background for context
        scatter(find(~isValid), history_Loss(~isValid), 20, [0.8 0.8 0.8], 'x'); 
        scatter(find(isValid), history_Loss(isValid), 30, [0.6 0.8 0.6], 'filled');
        
        % Plot the step line for the best valid loss found so far
        stairs(1:num_evals, best_valid_so_far, 'b-', 'LineWidth', 2.5);
        
        xlabel('Evaluation Number', 'FontWeight', 'bold');
        ylabel('Best Minimum Sats (Loss)', 'FontWeight', 'bold');
        title("Observed Best Valid Loss vs. Runs @ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid Runs', 'Valid Runs', 'Best Valid So Far', 'Location', 'northeast');
        grid on; hold off;
        
        exportgraphics(f3, fullfile(out_dir, 'Best_Loss_Over_Time.png'), 'Resolution', 300);
        close(f3);
        
        % --- PLOT 4: Architecture Map (Planes vs Sats per Plane) ---
        f4 = figure('Visible','off','Name', 'Architecture Map', 'Color', 'w'); hold on;
        
        % Jitter the points slightly since they are integers. 
        % Without jitter, an 8x8 failure and an 8x8 success would print directly on top of each other!
        jitter_x = planes + (rand(size(planes))-0.5)*0.4;
        jitter_y = sats_pp + (rand(size(sats_pp))-0.5)*0.4;
        
        % Plot invalid as faint red crosses
        scatter(jitter_x(~isValid), jitter_y(~isValid), 30, [0.8 0.3 0.3], 'x');
        % Plot valid as solid dots colored by how "cheap" they are
        scatter(jitter_x(isValid), jitter_y(isValid), 70, history_Loss(isValid), 'filled', 'MarkerEdgeColor', 'k');
        
        colormap('parula'); cb = colorbar; cb.Label.String = 'Num Sats';
        xlabel('Num Planes', 'FontWeight', 'bold');
        ylabel('Sats per Plane', 'FontWeight', 'bold');
        title("Valid Constellation Architectures@ " + num2str(Cfg.Orbit_height / 1000) + " km");
        legend('Invalid', 'Valid', 'Location', 'best');
        grid on; hold off;
        exportgraphics(f4, fullfile(out_dir, 'NumPlanes_SatsPerPlane.png'), 'Resolution', 300);
        close(f4);
        
        
        % --- PLOT 5: Orbital Mechanics (Phasing vs Inclination) ---
        % Does the phase shift only matter at certain inclinations?
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
        
        % --- PLOT 6: Parallel Coordinates (The "Top 5" Viewer) ---
        f6 = figure('Visible','off','Name', 'Parallel Coordinates (Top 5 Viewer)', 'Color', 'w');
        valid_data = history_X(isValid, :);
        valid_loss = history_Loss(isValid);
        
        if height(valid_data) > 0
            % 1. Inject the loss into the table
            valid_data.Total_Sats = valid_loss;
            
            % 2. Sort to find the Top 5 best runs
            [sorted_loss, sort_idx] = sort(valid_loss, 'ascend');
            num_to_highlight = min(5, length(sort_idx));
            
            % 3. Create a blank label list, defaulting to 'Other Valid Runs'
            category_labels = repmat({'Other Valid Runs'}, length(valid_loss), 1);
            cat_order = cell(num_to_highlight + 1, 1);
            
            % 4. Assign 'Rank 1', 'Rank 2', etc. to the best runs
            for i = 1:num_to_highlight
                label = sprintf('Rank %d (%d Sats)', i, sorted_loss(i));
                category_labels{sort_idx(i)} = label;
                cat_order{i} = label; % Keep track of the order for the legend
            end
            cat_order{end} = 'Other Valid Runs'; % Put 'Other' at the very bottom
            
            % 5. Add this new categorization to the table
            valid_data.Rank = categorical(category_labels, cat_order);
            
            % 6. Define EXACTLY which columns to draw as axes (This hides the 'Rank' column)
            coord_vars = {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees', 'Total_Sats'};
            
            % 7. Plot it, explicitly setting CoordinateVariables
            p = parallelplot(valid_data, 'CoordinateVariables', coord_vars, 'GroupVariable', 'Rank');
            
            % 8. Custom Colors: Bright distinct colors for Top 5, Light Gray for 'Other'
            custom_colors = lines(num_to_highlight); 
            custom_colors = [custom_colors; 0.85 0.85 0.85]; 
            p.Color = custom_colors;
            
            % 9. Custom Line Widths: An array with '4' for the Top 5, and '1' for 'Other'
            custom_widths = [repmat(4, 1, num_to_highlight), 1];
            p.LineWidth = custom_widths; 
            
            p.LineAlpha = 0.8; % Slight transparency so overlaps blend nicely
            
            title('Top 5 Solutions');
        else
            text(0.5, 0.5, 'No valid runs to plot.', 'HorizontalAlignment', 'center', 'FontSize', 14);
            axis off;
        end
        
        exportgraphics(f6, fullfile(out_dir, 'Top5_Solutions.png'), 'Resolution', 300);
        close(f6);
    
    
        % Show high res result of best constellation
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
        plot_results = true;
        use_parallel = true; 
        calc_link    = true;
        Cfg.StopTime   = datetime('2-Jun-2025 23:59:59', 'TimeZone', 'UTC'); % 48 hours
        Cfg.SampleTime = 20; % seconds
        
        Cfg.Lat_vec = linspace(55, 85, 10); % 32 workers on server
        Cfg.Lon_vec = linspace(-60, 30, 3);
        
        % Inject final parameters directly from surrogateopt's best_params
        Cfg.Num_planes     = best_params.Num_planes;
        Cfg.Sats_per_plane = best_params.Sats_per_plane;
        Cfg.Inclination    = best_params.Inclination;
        Cfg.Total_sats     = Cfg.Num_planes * Cfg.Sats_per_plane;
        
        % Calculate true Walker Phase constraint using modulo math
        Cfg.Phasing = mod(best_params.Phasing_Factor, Cfg.Num_planes);
        
        fprintf('\n--------------------------------------------------\n');
        fprintf('--- Simulating Optimal Constellation (Total Sats: %d) ---\n', Cfg.Total_sats);
        disp(best_params);
        
        % Run the detailed simulator!
        detailed_metrics = coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link);
        
        show_interactive = false;
        save_fig = true;
        show_constellation(Cfg, show_interactive, save_fig, out_dir);
     end
end

%% --- HELPER: Surrogate Optimization Wrapper ---
function obj_struct = surrogateWrapper(x_vec, Cfg)
    % 1. Inject variables into Cfg
    Cfg.Num_planes = x_vec(1);
    Cfg.Sats_per_plane = x_vec(2);
    Cfg.Inclination = x_vec(3);
    
    % 2. Calculate true Walker Phase constraint! (Integer modulo math)
    Cfg.Phasing = round((x_vec(3) / 360) * Cfg.Num_planes); 
    Cfg.Phasing = max(0, min(Cfg.Phasing, Cfg.Num_planes - 1)); 
    % Cfg.Phasing = mod(x_vec(4), Cfg.Num_planes); 
    
    Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
    
    % 3. Run your existing simulation
    metrics = coverage_simulator_function(Cfg, false, false, false);
    
    % 4. Package it into the struct surrogateopt demands
    obj_struct.Fval = Cfg.Total_sats;
    obj_struct.Ineq = 99.9 - metrics.worst_coverage_percent; 
end