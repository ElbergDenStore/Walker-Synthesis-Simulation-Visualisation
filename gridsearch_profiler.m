function [best_params, all_candidates] = gridsearch_profiler(Cfg, plot_results, min_sats)
    %% TURN ON PROFILER
    disp('Starting MATLAB Profiler...');
    profile on;
    
    %% 1. Build the Ascending Grid
    P_vec = 4:15; % Num Planes
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
    
    fprintf('\n=== Starting Ascending Grid Search (PROFILER MODE) ===\n');
    fprintf('Testing %d valid architectures from %d satellites...\n\n', height(search_grid), min_sats);
    
    %% 3. Simulate until we find the Global Minimum + N candidates
    best_params = [];
    all_candidates = table();
    
    % NOTE: Pool is NOT used. Profiler needs single-threaded execution to see inside coverage_simulator_function.
    num_workers = 1; 
    total_runs = height(search_grid);
    
    fprintf('Using %d SEQUENTIAL worker(s) to allow MATLAB Profiler to capture metrics...\n', num_workers);
    
    evaluated_coverage = NaN(total_runs, 1);
    status_flags = zeros(total_runs, 1); % 0 = Invalid, 1 = Candidate, 2 = Best
    
    min_sats_found = Inf;

    for batch_start = 1 : num_workers : total_runs
        if height(all_candidates) >= target_num_candidates
            fprintf('\n--- Found %d viable candidates. Target reached! Stopping search. ---\n', height(all_candidates));
            break;
        end
        
        batch_end = min(batch_start + num_workers - 1, total_runs);
        batch_size = batch_end - batch_start + 1;
        
        fprintf('\n--- Simulating Batch: Runs %d to %d ---\n', batch_start, batch_end);
        
        batch_coverage = zeros(batch_size, 1);
        batch_grid = search_grid(batch_start:batch_end, :);
        
        % --- THE SEQUENTIAL LOOP (CHANGED FROM PARFOR FOR PROFILING) ---
        for i = 1:batch_size
            local_Cfg = Cfg;
            local_Cfg.Num_planes     = batch_grid.Num_planes(i);
            local_Cfg.Sats_per_plane = batch_grid.Sats_per_plane(i);
            local_Cfg.Inclination    = batch_grid.Inclination(i);
            local_Cfg.Phasing        = batch_grid.Phasing(i); 
            local_Cfg.Total_sats     = batch_grid.Total_sats(i);
            
            % We will also force use_parallel=false so coverage_simulator_function runs sequentially
            metrics = coverage_simulator_function(local_Cfg, false, false, false);
            metrics = fast_coverage_simulator_function(local_Cfg, true, false, false); %reset cache to make sure new UES are correct
            metrics = fast_coverage_simulator_function(local_Cfg, false, false, false); %dont reset cache as UEs are exactly the same
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
            detailed_Cfg.Equal_UE_area = true; 
            
            detailed_Cfg.Num_planes     = temp_params.Num_planes;
            detailed_Cfg.Sats_per_plane = temp_params.Sats_per_plane;
            detailed_Cfg.Inclination    = temp_params.Inclination;
            detailed_Cfg.Phasing        = temp_params.Phasing;
            detailed_Cfg.Total_sats     = temp_params.Total_sats;
            
            fprintf('  -> High-Fidelity Test for %dx%d (Total: %d)...\n', ...
                detailed_Cfg.Num_planes, detailed_Cfg.Sats_per_plane, detailed_Cfg.Total_sats);
            
            detailed_metrics = coverage_simulator_function(detailed_Cfg, false, false, false); % plot_results=false, use_parallel=false
            
            
            if detailed_metrics.worst_coverage_percent > 99.999
                fprintf('  -> PASSED. Added to Candidates %d/%d.\n',height(all_candidates), target_num_candidates);
                
                status_flags(global_row_idx) = 1; 
                all_candidates = [all_candidates; temp_params];
                
                if temp_params.Total_sats < min_sats_found
                    min_sats_found = temp_params.Total_sats;
                    best_params = temp_params;
                    
                    fprintf('\n======================================================\n');
                    fprintf('New Minimum for %d km Found: %d Sats.\n', (Cfg.Orbit_height / 1000), min_sats_found);
                    fprintf('Current Candidates Pool: %d / %d\n', height(all_candidates), target_num_candidates);
                    fprintf('======================================================\n');
                end
            else
                fprintf('  -> Failed high-fidelity test (Cov: %.4f%%).\n', detailed_metrics.worst_coverage_percent);
            end
        end
        
        % Early exit capability for profile runs (stop after first 2 batches max)
        if batch_end >= 2
            fprintf('\n--- Early exit triggered in profiler to save time. Stopping after 2 runs ---\n');
            break;
        end
    end
    
    if isempty(best_params)
        fprintf('\n[!] PROFILER GRID SEARCH EXHAUSTED OR EARLY EXITED [!]\n');
    end
    
    %% --- END PROFILING & SAVE HTML ---
    fprintf('\nSaving MATLAB Profiler Results...\n');
    profile off;
    
    % Save to an HTML report folder named profile_results
    out_dir = fullfile('simulation_output', 'profiler_results');
    if exist(out_dir, 'dir')
        rmdir(out_dir, 's');
    end
    profsave(profile('info'), out_dir);
    
    fprintf('Profiling complete! Open `%s` in your browser or run `profile viewer` to inspect exactly where the time is being spent.\n', fullfile(out_dir, 'index.html'));
end
