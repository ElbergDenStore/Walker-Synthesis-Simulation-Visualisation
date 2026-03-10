% How to run through the night
% matlab -nodisplay -nosplash -nodesktop -batch "run_bayesian"
% -nodisplay -nosplash -nodesktop tells MATLAB to run purely as a command-line engine without booting up the heavy Java GUI interface.
% 
% -batch tells it to run your script, print the output directly to your terminal, and gracefully exit when it's done.
clear all; close all; clc;

%% 1. Define the Sweep Parameters
heights_km = 700:50:1200; % From 700 to 1200 in steps of 50
orbit_heights = heights_km * 1e3; % Convert to meters

num_runs = 100; % Number of iterations PER height 
plot_individual_results = true; % Set true to save the 6 detailed plots per height

% Preallocate an array to store the best satellite count for each height
% We use NaN (Not a Number) so we can easily skip heights that fail to find a valid solution
min_sats_array = NaN(size(orbit_heights)); 

dq = parallel.pool.DataQueue;
updateLiveScriptProgress(size(orbit_heights), true); 
afterEach(dq, @(~) updateLiveScriptProgress(size(orbit_heights), false));
start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
fprintf('\n=======================================================\n');
fprintf('SWEEP STARTED AT: %s\n', char(start_time));
fprintf('=======================================================\n');


%% 2. Run the Loop
for i = 1:length(orbit_heights)
    current_height = orbit_heights(i);
    
    fprintf('\n=======================================================\n');
    fprintf('Starting Optimization for Orbit Height: %d km (%d of %d)\n', ...
        heights_km(i), i, length(orbit_heights));
    fprintf('=======================================================\n');
    
    try
        % Call your optimizer function!
        best_params = bayesian_constellation_optimizer(current_height, num_runs, plot_individual_results);
        
        % Verify the optimizer actually returned a valid table row
        if ~isempty(best_params) && istable(best_params)
            % Calculate the total number of satellites from the winning parameters
            total_sats = best_params.Num_planes * best_params.Sats_per_plane;
            min_sats_array(i) = total_sats;
            fprintf('\n--> WINNER FOR %d km: %d Satellites\n', heights_km(i), total_sats);
        else
            fprintf('\n--> [!] No valid constellation found for %d km.\n', heights_km(i));
        end
        
    catch ME
        % If the optimizer crashes for one specific height, this prevents 
        % the entire multi-hour loop from terminating early!
        fprintf('\n[!] Error evaluating %d km: %s\n', heights_km(i), ME.message);
    end
    send(dq, []);
end

%% 3. Plot the Final Master Curve
figure('Name', 'Min Sats vs Orbit Height', 'Color', 'w');

% Only plot the heights where a valid constellation was successfully found
valid_idx = ~isnan(min_sats_array);

% Plot the curve with distinct markers
plot(heights_km(valid_idx), min_sats_array(valid_idx), '-ok', ...
    'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.7 0.2]);

xlabel('Orbit Altitude (km)', 'FontWeight', 'bold');
ylabel('Minimum Required Satellites', 'FontWeight', 'bold');
title('Optimal Constellation Size vs. Orbit Altitude (99.9% Coverage Requirement)');
grid on;

% Make the Y-axis strictly integers since you can't have half a satellite
yticks(min(min_sats_array(valid_idx)):max(min_sats_array(valid_idx)));

% Save the master plot and the raw data array so you don't lose the results!
save('Master_Altitude_Sweep_Results.mat', 'heights_km', 'min_sats_array');
exportgraphics(gcf, 'Master_MinSats_vs_Altitude.png', 'Resolution', 300);

fprintf('\n=== SWEEP COMPLETE ===\n');
fprintf('Master plot saved as Master_MinSats_vs_Altitude.png\n');

% --- RECORD END TIME ---
end_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
elapsed_time = end_time - start_time;

fprintf('\n=======================================================\n');
fprintf('SWEEP FINISHED AT: %s\n', char(end_time));
fprintf('TOTAL ELAPSED TIME: %s\n', char(elapsed_time));
fprintf('Master plot saved as Master_MinSats_vs_Altitude.png\n');
fprintf('=======================================================\n');


function updateLiveScriptProgress(total_pts, reset_flag)
    persistent p last_percent reverseStr
    
    % Initialization / Reset
    if nargin > 1 && reset_flag
        p = 0;
        last_percent = -1; 
        reverseStr = '';
        return;
    end
    
    if isempty(p)
        p = 0;
        last_percent = -1;
        reverseStr = '';
    end
    
    p = p + 1;
    current_percent = floor((p / total_pts) * 100);
    
    % Only update the screen when the percentage actually changes (prevents terminal lag)
    if current_percent > last_percent || p == total_pts
        bar_length = 40; % How wide you want the progress bar to be
        num_equals = round((current_percent / 100) * bar_length);
        num_spaces = bar_length - num_equals;
        
        % Build the string: e.g., [========          ]
        bar_str = ['[', repmat('=', 1, num_equals), repmat(' ', 1, num_spaces), ']'];
        
        % Create the full message
        msg = sprintf('Processing: %s %d%%', bar_str, current_percent);
        
        % Print backspaces to clear the old line, then print the new line
        fprintf([reverseStr, msg]);
        
        % Save the number of backspaces needed for the next loop
        reverseStr = repmat(sprintf('\b'), 1, length(msg)-1);
        
        last_percent = current_percent;
    end
    
    % Cap it off cleanly when finished and drop to a new line
    if p >= total_pts
        % fprintf('\n');
        p = 0; 
        last_percent = -1;
        reverseStr = ''; 
    end
end