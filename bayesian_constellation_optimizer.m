function best_params = bayesian_constellation_optimizer(orbit_height, num_runs, plot_results)

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
    
    search_vars = [var_planes, var_sats, var_inc, var_phasing];
    
    %% Bayesian Optimization
    obj_fun = @(x) evaluateConstellationLoss(x, Cfg); 
    
    % Run optimizer
    results = bayesopt(obj_fun, search_vars, ...
        'MaxObjectiveEvaluations', num_runs, ... 
        'IsObjectiveDeterministic', true, ...
        'ExplorationRatio', 0.5, ...
        'UseParallel', true, ...
        'NumCoupledConstraints', 1, ... 
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'PlotFcn', {@plotObjectiveModel, @plotMinObjective});
    
    %% Extract and Save the Best Results
    best_params = results.XAtMinObjective;
    fprintf('\n=== OPTIMIZATION COMPLETE ===\n');
    fprintf('Best Constellation Found:\n');
    disp(best_params);
    
    % Pull the actual metrics for the BEST run, not the first run
    custom_metrics_history = results.UserDataTrace;
    best_idx = results.IndexOfMinimumTrace(end); 
    
    fprintf('--- Metrics for the Best Valid Constellation ---\n');
    disp(custom_metrics_history{best_idx});
    
    
    if plot_results
        %% Create Output Directory
        % Format: simulation_output/bayesian_runs/orbit_sats_inclination_phasing_date
        date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name = sprintf('%.0f_%d', ...
            Cfg.Orbit_height/1e3, custom_metrics_history{best_idx}.Total_sats, date_str);
        out_dir = fullfile('simulation_output/bayesian_runs', folder_name);
        
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end

        %% --- POST-RUN VISUALIZATIONS ---
        history_X = results.XTrace; 
        history_Loss = results.ObjectiveTrace;
        history_Constraints = results.ConstraintsTrace;
        
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
        cov_history = cellfun(@(x) x.Cov_percent, results.UserDataTrace);
        sat_history = cellfun(@(x) x.Total_sats, results.UserDataTrace);
        
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
    end 
    
    % Show high res result of best constellation

    valid_indices = find(results.ConstraintsTrace <= 0);

    if isempty(valid_indices)
        fprintf('[!] No valid constellations found to simulate.\n');
    else
        % 2. Get the Total Satellites (Objective) for those valid runs
        valid_objectives = results.ObjectiveTrace(valid_indices);
        
        % 3. Sort them from lowest number of sats to highest
        [sorted_sats, sort_order] = sort(valid_objectives, 'ascend');
        
        % 4. Map the sorted list back to the original run numbers
        best_valid_indices = valid_indices(sort_order);
        
        % 5. Cap amount of runs
        num_to_run = min(1, length(best_valid_indices));
        fprintf('Found %d valid constellations. Running detailed simulations for the top %d...\n\n', length(best_valid_indices), num_to_run);

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
        
        for i = 1:num_to_run
            % Extract the parameters for this specific rank
            idx = best_valid_indices(i);
            best_params = results.XTrace(idx, :);
            
            fprintf('\n--------------------------------------------------\n');
            fprintf('--- Simulating Rank %d (Total Sats: %d, Run Index: %d) ---\n', i, sorted_sats(i), idx);
            disp(best_params);
            
            % 7. Inject these parameters back into Cfg
            Cfg.Num_planes     = best_params.Num_planes;
            Cfg.Sats_per_plane = best_params.Sats_per_plane;
            Cfg.Inclination    = best_params.Inclination;
            
            % Do the exact same Phasing math we did in the objective function
            Cfg.Phasing = round((best_params.Phasing_Degrees / 360) * Cfg.Num_planes);
            Cfg.Phasing = max(0, min(Cfg.Phasing, Cfg.Num_planes - 1));
            Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
            
            % 8. Run the simulator!
            detailed_metrics = coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link);
            
        end
    end
end


%% --- HELPER: The Objective Function ---
function [loss, constraints, UserData] = evaluateConstellationLoss(params, Cfg)
    opt_vars = params.Properties.VariableNames;
    
    for k = 1:length(opt_vars)
        var_name = opt_vars{k};
        Cfg.(var_name) = params.(var_name);
    end
    
    Cfg.Phasing = round((params.Phasing_Degrees / 360) * Cfg.Num_planes); 
    Cfg.Phasing = max(0, min(Cfg.Phasing, Cfg.Num_planes - 1)); 
    Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
    
    plot_results = false;
    use_parallel = false;
    calc_link = false;
    
    % Get metrics from simulator
    metrics = coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link);
    
    % 1. Calculate Loss (Minimize total satellites)
    loss = Cfg.Total_sats;
    
    % 2. Constraints (Require >= 99.9% coverage)
    constraints = 99.9 - metrics.worst_coverage_percent; 

    % 3. Package metrics for later analysis
    UserData.Cov_percent  = metrics.worst_coverage_percent;
    UserData.Total_sats   = Cfg.Total_sats; 
    % UserData.Num_visible  = metrics.Num_visible; 
end