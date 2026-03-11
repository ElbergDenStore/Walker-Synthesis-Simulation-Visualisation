% How to run through the night
% matlab -nodisplay -nosplash -nodesktop -batch "run_bayesian"
% -nodisplay -nosplash -nodesktop tells MATLAB to run purely as a command-line engine without booting up the heavy Java GUI interface.
% 
% -batch tells it to run your script, print the output directly to your terminal, and gracefully exit when it's done.
% Press Ctrl+b, release both keys, and then press d.
clear all; close all; clc;

%% 1. Define the Sweep Parameters
heights_km = 700:25:1200; % From 700 to 1200 in steps of 50
orbit_heights = heights_km * 1e3; % Convert to meters

num_runs = 1000; % Number of iterations PER height 
plot_individual_results = false; % Set true to save the 6 detailed plots per height

% Preallocate an array to store the best satellite count for each height
% We use NaN (Not a Number) so we can easily skip heights that fail to find a valid solution
min_sats_array = NaN(size(orbit_heights)); 

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
        % best_params = bayesian_constellation_optimizer(current_height, num_runs, plot_individual_results);
        best_params = surrogateopt_constellation_optimizer(current_height, num_runs, plot_individual_results);
        % best_params = gridsearch(current_height);

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
end

%% 3. Plot the Final Master Curvefigure('Name', 'Min Sats vs Orbit Height', 'Color', 'w');

% Only plot the heights where a valid constellation was successfully found
valid_idx = ~isnan(min_sats_array);
f1 = figure('Visible', 'off', 'Name', 'Bayesian_sweep_results', 'Color', 'w', 'Position', [100 100 1000 600]); 
% Plot the curve with distinct markers
scatter(heights_km(valid_idx), min_sats_array(valid_idx),25, 'filled', ...
                'MarkerFaceColor', '#0072BD', 'MarkerFaceAlpha', 1);

xlabel('Orbit Height (km)', 'FontWeight', 'bold');
ylabel('Num Sats', 'FontWeight', 'bold');
title('Optimal Constellation and Orbit Height');
grid on;

% Make the Y-axis strictly integers since you can't have half a satellite
% yticks(min(min_sats_array(valid_idx)):max(min_sats_array(valid_idx)));
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
folder_name = sprintf('Bayesian_sweep_%s', date_str);
out_dir = fullfile('simulation_output', folder_name);
            
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Save the master plot and the raw data array so you don't lose the results!
save(fullfile(out_dir,'Master_Altitude_Sweep_Results.mat'),'heights_km', 'min_sats_array');
exportgraphics(f1, fullfile(out_dir, 'Bayesian_sweep_results.png'), 'Resolution', 300);
close(f1);
% exportgraphics(gcf, 'Master_MinSats_vs_Altitude.png', 'Resolution', 300);

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
