function plot_sweep(sweep_folder)
% PLOT_SWEEP  Regenerate master sweep plot from saved results.
%
% Usage:
%   plot_sweep()                   - most recent Master_Sweep_* folder
%   plot_sweep('path/to/sweep')    - specific sweep folder

script_dir     = fileparts(mfilename('fullpath'));
workspace_root = fileparts(script_dir);
sim_dir        = fullfile(workspace_root, 'simulation_output');

if nargin == 0
    hits = dir(fullfile(sim_dir, 'Master_Sweep_*'));
    hits = hits([hits.isdir]);
    if isempty(hits)
        error('No Master_Sweep_* folders found in:\n  %s', sim_dir);
    end
    [~, idx]     = max([hits.datenum]);
    sweep_folder = fullfile(sim_dir, hits(idx).name);
    fprintf('Using most recent sweep: %s\n', sweep_folder);
else
    % Resolve relative paths against workspace root
    sweep_folder = char(sweep_folder);
    if sweep_folder(1) ~= '/' && ~(numel(sweep_folder) > 1 && sweep_folder(2) == ':')
        sweep_folder = fullfile(workspace_root, sweep_folder);
    end
end

mat_file = fullfile(sweep_folder, 'Master_Altitude_Sweep_Results.mat');
if ~isfile(mat_file)
    error('Master_Altitude_Sweep_Results.mat not found in:\n  %s', sweep_folder);
end

loaded = load(mat_file);

heights_km      = loaded.heights_km;
star_sats       = loaded.star_sats;
best_delta_sats = loaded.best_delta_sats;

star_plot_y  = [star_sats.Total_sats];
delta_plot_y = best_delta_sats.Total_sats';

if isfield(loaded, 'Master_config')
    minimum_lat_deg = min(loaded.Master_config.Lat_range_deg);
    title_str = sprintf('Walker Star vs Minimum Walker Delta (Lat: %0.1f\x00B0)', minimum_lat_deg);
else
    title_str = 'Walker Star vs Minimum Walker Delta';
end

set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

f1 = figure('Visible', 'off', 'Name', 'Constellation Comparison', 'Color', 'w', 'Position', [100 100 1000 600]);
hold on;
scatter(heights_km, star_plot_y,  36, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'Analytical Walker Star');
scatter(heights_km, delta_plot_y, 36, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Optimized Walker Delta');
xlabel('Orbit Height (km)', 'FontWeight', 'bold');
ylabel('Total Satellites Required', 'FontWeight', 'bold');
title(title_str);
legend('Location', 'northeast');
grid on; hold off;

out_path = fullfile(sweep_folder, 'Star_vs_Delta_Comparison.png');
exportgraphics(f1, out_path, 'Resolution', 300);
close(f1);
fprintf('Plot saved to: %s\n', out_path);
end
