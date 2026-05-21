% Walker Star gridsearch sweep – finds the minimum-satellite Walker Star
% constellation that achieves 99.999% coverage at each altitude.
%
% Inclination is fixed at 87°, phasing at floor(P/2).
% Search space: 2–10 planes, 5–30 sats/plane.
%
% How to run through the night:
% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "Find_valid_walker_star_gridsearch" > log.txt

clear; close all; clc;
delete(gcp('nocreate'));

%% Master Configuration
heights_km                          = 500:10:1200;
Master_config.Lat_range_deg         = [54+(35/60), 83+(40/60)];
Master_config.Min_elevation_UE      = 20;
Master_config.Num_Planes            = 2:15;
Master_config.Sats_Plane            = 5:35;
Master_config.Inc_vec               = 90; 
Master_config.WalkerStar            = true; 
Master_config.Target_num_candidates = 1;
Master_config.SampleTime            = 660;

% Sub Run configurations (same statistical guarantees as delta sweep)
Master_config.Ultrafast.Duration_h  = 3;
Master_config.Ultrafast.Num_UEs     = 200;
Master_config.Fast.Duration_h       = 50;
Master_config.Fast.Num_UEs          = 800;
certainty       = 99   * 1e-2;
fractional_area = 0.1  * 1e-2;
fractional_time = 0.1  * 1e-2;
required_samples = log(1-certainty) / log(1 - fractional_area*fractional_time)

required_time_h = ceil(sqrt(required_samples)) / (3600/Master_config.SampleTime)
Master_config.Detailed.Duration_h  = required_time_h;
required_UEs = ceil(sqrt(required_samples))
Master_config.Detailed.Num_UEs     = required_UEs;

% Record start time
start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
fprintf('=======================================================\n');
fprintf('STARTING WALKER STAR GRIDSEARCH SWEEP AT: %s\n', char(start_time));
fprintf('=======================================================\n');

% Create master sweep output folder up-front so gridsearch runs nest inside it
date_str_start  = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
folder_name     = sprintf('Master_Sweep_WalkerStar_%s', date_str_start);
out_dir         = fullfile('simulation_output', folder_name);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
gridsearch_base = fullfile(out_dir, 'gridsearch_runs');

star_sats       = [];
best_star_sats  = table();
all_star_sats   = {};
gridsearch_dirs = cell(length(heights_km), 1);

% Descending order so min_sats hint propagates from higher (easier) altitudes
heights_km = sort(heights_km, 'descend');
for i = 1:length(heights_km)

    %% Analytical Walker Star as reference
    minimum_lat_deg = min(Master_config.Lat_range_deg);
    [Num_planes, Sats_per_plane, Total_sats] = calculate_walker_star( ...
        heights_km(i), minimum_lat_deg, Master_config.Min_elevation_UE);
    star_sats(i).Orbit_height   = heights_km(i);
    star_sats(i).Total_sats     = Total_sats;
    star_sats(i).Num_planes     = Num_planes;
    star_sats(i).Phasing        = Num_planes/2;
    star_sats(i).Inclination    = 90;
    star_sats(i).Sats_per_plane = Sats_per_plane;

    if i > 1
        min_sats = best_star_sats.Total_sats(i-1);
    else
        min_sats = 0;
    end

    [best_params, all_star_sats{i}, gridsearch_dirs{i}] = gridsearch( ...
        Master_config, heights_km(i), min_sats, gridsearch_base);

    best_star_sats = [best_star_sats; best_params(1, :)];
end

%% Save Outputs
save(fullfile(out_dir, 'Master_Altitude_Sweep_Results.mat'), ...
    'heights_km', 'star_sats', 'best_star_sats', 'all_star_sats', ...
    'Master_config', 'gridsearch_dirs');
fprintf('Results saved to: %s\n', out_dir);

end_time     = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
elapsed_time = end_time - start_time;
fprintf('\n=======================================================\n');
fprintf('SWEEP FINISHED AT: %s\n', char(end_time));
fprintf('TOTAL ELAPSED TIME: %s\n', char(elapsed_time));
fprintf('=======================================================\n');
