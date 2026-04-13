% How to run through the night:
% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "Find_optimal_constellations" > log.txt
% xvfb is a virtual display to avoid constellation pictures do not crash server
% one ">" overwrites the file
% Read log during run using tail -f log.txt

% Read log afterwards using less
% How to use it: Type less log.txt
% Pro-tips inside less:% 
% Press Space to page down, b to page up.
% Press G to jump immediately to the very bottom (the newest logs).
% Press g to jump back to the top.
% Type / followed by a keyword (like /error or /crash) and hit Enter to search. Press n to jump to the next match
% Press q to quit.

clear; close all; clc;
%% 1. Force kill the current parallel pool
    delete(gcp('nocreate')); % necessary or it will get stuck

%% Master Configuration
heights_km                          = 700:10:1200;
plot_individual_results = true;
Master_config.Lat_range_deg         = [55, 85];
Master_config.Min_elevation_UE      = 20;
Master_config.Num_Planes            = 2:20; % Num Planes
Master_config.Sats_Plane            = 2:20; % Sats per Plane
Master_config.Inc_vec               = linspace(70, 80, 21); % More general -> linspace(max(Lat_range_deg)-15, min(max(Lat_range_deg),80), 21)
Master_config.Target_num_candidates = 10;

% Sub Run configurations
Master_config.Ultrafast.Duration_h  = 1;  
Master_config.Ultrafast.Num_UEs     = 100;
Master_config.Fast.Duration_h       = 24;  
Master_config.Fast.Num_UEs          = 100;
Master_config.Detailed.Duration_h   = 36;      
Master_config.Detailed.Num_UEs      = 2000;


% Record start time for the sweep
start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
fprintf('=======================================================\n');
fprintf('STARTING CONSTELLATION SWEEP AT: %s\n', char(start_time));
fprintf('=======================================================\n');


star_sats = [];
best_delta_sats = table();
all_delta_sats = [];

% Runs through all orbit heights from top to bottom
heights_km = sort(heights_km,"descend");
for i = 1:length(heights_km) 
    current_h_meters = heights_km(i) * 1000;
  
    %% The Analytical "Seed" (Walker Star Baseline)
    [Num_planes, Sats_per_plane, Total_sats] = calculate_walker_star(heights_km(i), min(Master_config.Lat_range_deg), Master_config.Min_elevation_UE);
    star_sats(i).Orbit_height   = heights_km(i);
    star_sats(i).Total_sats     = Total_sats;
    star_sats(i).Num_planes     = Num_planes;
    star_sats(i).Phasing        = Num_planes/2;
    star_sats(i).Inclination    = 87;
    star_sats(i).Sats_per_plane = Sats_per_plane;


    if i > 1
        min_sats = best_delta_sats.Total_sats(i-1); % Limit search space based on previous result
    else
        min_sats = 0; %floor(star_N * 0.6); % Limit search space based on analytical star
    end 

    [best_params, all_delta_sats{i}] = gridsearch(Master_config, heights_km(i), plot_individual_results, min_sats);

    % Keep one best row per altitude in sweep order.
    best_delta_sats = [best_delta_sats; best_params(1, :)];
    %% Save the Optimized Result
    % if ~isempty(best_params)
    %     best_delta_sats(i).Orbit_height = heights_km(i);
    %     best_delta_sats(i).Total_sats     = best_params.Total_sats;
    %     best_delta_sats(i).Num_planes     = best_params.Num_planes;
    %     best_delta_sats(i).Phasing        = best_params.Phasing;
    %     best_delta_sats(i).Sats_per_plane = best_params.Sats_per_plane;
    %     best_delta_sats(i).Inclination    = best_params.Inclination;
    % else
    %     % if no solutions are found, it deserves to crash
    %     error_msg = sprintf('FATAL ERROR: Optimizer failed to find a solution at %d km! Halting sweep.', heights_km(i));
    %     error(error_msg);
    % end
end

%% --- 3. Plot the Final Master Curve (Star vs Delta) ---
star_plot_y  = [star_sats.Total_sats];
delta_plot_y = [best_delta_sats.Total_sats];
f1 = figure('Visible', 'off', 'Name', 'Constellation Comparison', 'Color', 'w', 'Position', [100 100 1000 600]); hold on;

% Plot the Analytical Walker Star baseline (Red Line)
plot(heights_km, star_plot_y, '-ro', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Analytical Walker Star');

% Plot the Optimized Walker Delta results (Blue Line)
% We only plot valid indices in case one of the heights failed
plot(heights_km, delta_plot_y, '-bs', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Optimized Walker Delta');

xlabel('Orbit Height (km)', 'FontWeight', 'bold');
ylabel('Total Satellites Required', 'FontWeight', 'bold');
title(sprintf('Coverage Efficiency: Walker Star vs Walker Delta (Lat: %d°)', target_lat));
legend('Location', 'northeast');
grid on; hold off;

%% --- 4. Save Outputs ---
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
folder_name = sprintf('Master_Sweep_%s', date_str);
out_dir = fullfile('simulation_output', folder_name);
            
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Save the data and the plot
save(fullfile(out_dir,'Master_Altitude_Sweep_Results.mat'), 'heights_km', 'star_sats', 'best_delta_sats','all_delta_sats');
exportgraphics(f1, fullfile(out_dir, 'Star_vs_Delta_Comparison.png'), 'Resolution', 300);
close(f1);

% --- RECORD END TIME ---
end_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
elapsed_time = end_time - start_time;

fprintf('\n=======================================================\n');
fprintf('SWEEP FINISHED AT: %s\n', char(end_time));
fprintf('TOTAL ELAPSED TIME: %s\n', char(elapsed_time));
fprintf('Master plot saved as Star_vs_Delta_Comparison.png\n');
fprintf('=======================================================\n');
