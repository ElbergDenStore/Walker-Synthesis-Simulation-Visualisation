% How to run through the night
% matlab -nodisplay -nosplash -nodesktop -batch "run_bayesian"
% -nodisplay -nosplash -nodesktop tells MATLAB to run purely as a command-line engine without booting up the heavy Java GUI interface.
% 
% -batch tells it to run your script, print the output directly to your terminal, and gracefully exit when it's done.
% Press Ctrl+b, release both keys, and then press d.

% function best_params = master_constellation_optimizer(orbit_height,
% method, plot_results) - might make it into a function

clear all; close all; clc;
method = "grid"
orbit_height = 700e3; % can make it a vector
target_lat = 55;
plot_results = true;


%% 1. Base Configuration
Cfg.StartTime  = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
Cfg.StopTime   = datetime('1-Jun-2025 11:59:59', 'TimeZone', 'UTC'); 
Cfg.SampleTime = 60; 
Cfg.Lat_vec = linspace(55, 85, 6); 
Cfg.Lon_vec = linspace(-60, 30, 1);
Cfg.Min_elevation_UE = 20;
Cfg.WalkerStar     = false; % optimizing walker deltas
Cfg.Orbit_height = orbit_height;

%% 2. The Analytical "Seed" (Walker Star Baseline)
% Call your custom math to find the theoretical upper bound
[star_P, star_S, star_N] = get_analytical_star(orbit_height / 1000, target_lat, Cfg.Min_elevation_UE);

fprintf('\n======================================================\n');
fprintf('Analytical Star Baseline: %d Planes x %d Sats (%d Total)\n', star_P, star_S, star_N);
fprintf('======================================================\n');

%% 3. Dynamically Bound the Search Space
min_sats = floor(star_N * 0.7); % Search down to 70%
max_sats = star_N*1.2;          % max


%% 4. Route to the chosen Optimizer
switch lower(method)
    case 'grid'
        fprintf('Running Smart Ascending Grid Search...\n');
        best_params = gridsearch(Cfg,plot_results, min_sats, max_sats);
        
    case 'surrogate'
        fprintf('Running Surrogate Optimization...\n');
        best_params = run_surrogate(Cfg, num_runs, min_sats, max_sats);
        
    case 'bayes'
        fprintf('Running Bayesian Optimization...\n');
        best_params = run_bayes(Cfg, num_runs, min_sats, max_sats);
        
    otherwise
        error('Invalid method. Choose: ''grid'', ''surrogate'', or ''bayes''.');
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