% run_altitude_sweep.m
clear all; close all; clc;

%% 1. Define the Sweep Parameters
heights_km = 700:50:1200; % From 700 to 1200 in steps of 50
orbit_heights = heights_km * 1e3; % Convert to meters

num_runs = 100; % Number of iterations PER height 
plot_individual_results = true; % Set true to save the 6 detailed plots per height

% Preallocate an array to store the best satellite count for each height
% We use NaN (Not a Number) so we can easily skip heights that fail to find a valid solution
min_sats_array = NaN(size(orbit_heights)); 

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