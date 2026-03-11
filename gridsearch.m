function best_params = gridsearch(orbit_height)
    %% 1. Build the Ascending Grid
    P_vec = 4:15; % Num Planes
    S_vec = 4:15; % Sats per Plane
    Inc_vec = linspace(70, 80, 11); % 6 Inclinations
    
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
    % Filter out anything that isn't between 50 and 100 satellites
    isValidTarget = search_grid.Total_Sats >= 50 & search_grid.Total_Sats <= 100;
    search_grid = search_grid(isValidTarget, :);
    
    % Sort from cheapest to most expensive!
    search_grid = sortrows(search_grid, 'Total_Sats', 'ascend');
    
    fprintf('\n=== Starting Ascending Grid Search ===\n');
    fprintf('Testing %d valid architectures between 50 and 100 satellites...\n\n', height(search_grid));
    
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
            % Initialize Cfg completely inside the loop for parallel safety
            Cfg = struct();
            Cfg.StartTime  = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
            Cfg.StopTime   = datetime('1-Jun-2025 11:59:59', 'TimeZone', 'UTC'); 
            Cfg.SampleTime = 60; 
            Cfg.Lat_vec = linspace(55, 85, 6); 
            Cfg.Lon_vec = linspace(-60, 30, 1);
            Cfg.Min_elevation_UE = 20;
            Cfg.WalkerStar     = false;
            Cfg.Orbit_height   = orbit_height;
            
            % Inject architecture parameters
            Cfg.Num_planes     = batch_grid.Num_planes(i);
            Cfg.Sats_per_plane = batch_grid.Sats_per_plane(i);
            Cfg.Inclination    = batch_grid.Inclination(i);
            Cfg.Phasing        = batch_grid.Phasing_Factor(i); 
            Cfg.Total_sats     = batch_grid.Total_Sats(i);
            
            % Run the simulator!
            metrics = coverage_simulator_function(Cfg, false, false, false);
            batch_coverage(i) = metrics.worst_coverage_percent;
            
            fprintf('Finished %dx%d (Inc: %.1f, Phase: %d) -> Cov: %.2f%%\n', ...
                Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Inclination, Cfg.Phasing, batch_coverage(i));
        end
        % --- END PARALLEL LOOP ---
        
        % Evaluate the batch results
        valid_indices = find(batch_coverage >= 99.9);
        
        if ~isempty(valid_indices)
            % Because the grid is pre-sorted by satellite count, 
            % the FIRST index in valid_indices is guaranteed to be the cheapest!
            best_local_idx = valid_indices(1);
            best_params = batch_grid(best_local_idx, :);
            
            fprintf('\n======================================\n');
            fprintf('=== EUREKA! GLOBAL MINIMUM FOUND ===\n');
            fprintf('======================================\n');
            fprintf('Total Satellites: %d\n', best_params.Total_Sats);
            disp(best_params);
            
            % BREAK THE OUTER LOOP! We are done.
            break; 
        end
    end
end