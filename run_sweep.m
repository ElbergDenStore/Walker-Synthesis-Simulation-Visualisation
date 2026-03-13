% How to run through the night:
% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run_sweep" > sweep_log.txt
% xvfb is a virtual display to avoid constellation pictures do not crash
% server
% one ">" overwrites the file
% Read log during run using tail -f sweep_log.txt

% Read log afterwards using less
% How to use it: Type less my_log_file.txt.
% Pro-tips inside less:% 
% Press Space to page down, b to page up.
% Press G to jump immediately to the very bottom (the newest logs).
% Press g to jump back to the top.
% Type / followed by a keyword (like /error or /crash) and hit Enter to search. Press n to jump to the next match
% Press q to quit.

clear; close all; clc;

%% --- 1. Master Configuration ---
method = "grid"; % Toggle: 'grid', 'surrogate', or 'bayes'
heights_km = 700:10:1200; % Iterate over these altitudes (in km)
target_lat = 55;
plot_individual_results = true; % Keep false for the sweep to save time

% Record start time for the sweep
start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
fprintf('=======================================================\n');
fprintf('STARTING CONSTELLATION SWEEP AT: %s\n', char(start_time));
fprintf('=======================================================\n');

% Arrays to store the data for our combined plot
% star_sats_array = NaN(size(heights_km));
% delta_sats_array = NaN(size(heights_km));

star_sats = [];
best_delta_sats = [];
all_delta_sats = [];
%% --- 2. The Sweep Loop ---
for i = 1:length(heights_km)
    current_h_meters = heights_km(i) * 1000;
    
    %% Base Cfg for this iteration
    Cfg.StartTime  = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
    Cfg.StopTime   = datetime('1-Jun-2025 01:59:59', 'TimeZone', 'UTC'); 
    Cfg.SampleTime = 60; 
    Cfg.Lat_vec = linspace(55, 85, 4); 
    Cfg.Lon_vec = linspace(-60, 30, 1);
    Cfg.Min_elevation_UE = 20;
    Cfg.WalkerStar     = false; % We are optimizing Walker Deltas
    Cfg.Orbit_height = current_h_meters;

    %% The Analytical "Seed" (Walker Star Baseline)
    [star_P, star_S, star_N] = get_analytical_star(heights_km(i), target_lat, Cfg.Min_elevation_UE);
    star_sats(i).Orbit_height = heights_km(i);
    star_sats(i).Num_sats = star_N;
    star_sats(i).Num_planes = star_P;
    star_sats(i).Phasing = star_P/2;
    star_sats(i).Inclination = 87;
    star_sats(i).Sats_per_plane = star_S;

    fprintf('\n======================================================\n');
    fprintf('ALTITUDE: %d km\n', heights_km(i));
    fprintf('Analytical Star Baseline: %d Planes x %d Sats (%d Total)\n', star_P, star_S, star_N);
    fprintf('======================================================\n');

    %% Dynamically Bound the Search Space
    % We expect the Delta to beat the Star, so we look between 70% and 120% of the Star's sats
    min_sats = floor(star_N * 0.6); 
    max_sats = ceil(star_N);          

    %% Route to the chosen Optimizer
    best_params = [];
    switch lower(method)
        case 'grid'
            fprintf('Running Smart Ascending Grid Search...\n');
            % NOTE: Ensure your gridsearch function accepts these inputs!
            [best_params, all_delta_sats{i}] = gridsearch(Cfg, plot_individual_results, min_sats, max_sats);
            
        case 'surrogate'
            fprintf('Running Surrogate Optimization...\n');
            best_params = run_surrogate(Cfg, plot_individual_results, min_sats, max_sats);
            
        case 'bayes'
            fprintf('Running Bayesian Optimization...\n');
            best_params = run_bayes(Cfg, plot_individual_results, min_sats, max_sats);
            
        otherwise
            error('Invalid method. Choose: ''grid'', ''surrogate'', or ''bayes''.');
    end
    
    %% Save the Optimized Result
    if ~isempty(best_params)
        best_delta_sats(i).Orbit_height = heights_km(i);
        best_delta_sats(i).Num_sats       = best_params.Total_Sats;
        best_delta_sats(i).Num_planes     = best_params.Num_planes;
        best_delta_sats(i).Phasing_Factor = best_params.Phasing_Factor;
        best_delta_sats(i).Sats_per_plane = best_params.Sats_per_plane;
        best_delta_sats(i).Inclination    = best_params.Inclination;
    else
        % if no solutions are found, it deserves to crash
        error_msg = sprintf('FATAL ERROR: Optimizer failed to find a solution at %d km! Halting sweep.', heights_km(i));
        error(error_msg);
    end
end

%% --- 3. Plot the Final Master Curve (Star vs Delta) ---
star_plot_y  = [star_sats.Num_sats];
delta_plot_y = [best_delta_sats.Num_sats];
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


%% --- HELPER FUNCTIONS ---
function [P, S, N] = get_analytical_star(alt_km, phi_min, eps_min)
    Re = 6378.137;           
    Rs = Re + alt_km; 
    earth_O_at_lat = (cosd(phi_min)*Re)*2*pi;
    
    alpha = asind((Re / Rs) * cosd(eps_min));
    ECA = deg2rad(180 - (90 + eps_min + alpha));
    
    hex_side = ECA * Re;
    P = ceil((((earth_O_at_lat / 2) - hex_side) / (1.5 * hex_side)) + 1);
    S = ceil((2*pi)/(sqrt(3)*ECA)); 
    N = P * S;
end