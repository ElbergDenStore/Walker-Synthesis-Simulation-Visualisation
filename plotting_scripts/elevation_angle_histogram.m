clearvars; close all; clc;

%% 1) Simulation configuration
height_km = 1000;
constellation_type = "walkerdelta";
ue_grid_size = "big";
duration = "long";  % use "long" for a larger/full run
frequency = "ku";

use_parallel = false;
calc_link = false;     % Geometry only is enough for elevation histogram

Cfg = get_cfg(height_km, constellation_type, ue_grid_size, duration, frequency);

%% 2) Run simulation
fprintf('Running simulation for elevation histogram...\n');
metrics = coverage_simulator_function(Cfg, use_parallel, calc_link);

%% 3) Collect elevation samples across all UEs and time steps
SimDataArray = [metrics.SimData];
all_el_deg = vertcat(SimDataArray.Elevation_deg);
all_el_deg = all_el_deg(isfinite(all_el_deg));

if isempty(all_el_deg)
	warning('No valid elevation samples found. Nothing to plot.');
	return;
end

%% 4) Plot histogram
f = figure('Color', 'w', 'Position', [100 100 1000 550]);
histogram(all_el_deg, 50, 'Normalization', 'pdf', 'FaceColor', [0.0 0.45 0.74], 'EdgeColor', 'none');
grid on; box on;

xlabel('Elevation Angle (deg)', 'FontWeight', 'bold');
ylabel('Probability Density', 'FontWeight', 'bold');
title(sprintf('Elevation Angle Distribution | %d km, min el = %d', ...
	round(Cfg.Orbit_height/1e3),Cfg.Min_elevation_UE), 'FontSize', 14);

%% 5) Save at 600 DPI
out_dir = fullfile('figures', 'elevation_angle_histogram');
if ~exist(out_dir, 'dir')
	mkdir(out_dir);
end

% timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_file = fullfile(out_dir, sprintf('elevation_hist.png'));
exportgraphics(f, out_file, 'Resolution', 600);

fprintf('Saved figure: %s\n', out_file);
