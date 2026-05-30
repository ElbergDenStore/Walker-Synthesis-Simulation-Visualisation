clear; close all; clc;

% Ensure functions/ subdirectory is on the MATLAB path
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

%% Parameters
heights_km       = 500:0.1:1200;
Min_latitude_deg = 54 + 35/60;   % 54°35'N — southernmost Denmark
Min_elevation_UE = 20;            % degrees minimum elevation angle
inclination_deg  = 90;            % inclination to evaluate

n_heights = numel(heights_km);

%% Compute results
total_sats     = zeros(1, n_heights);
planes         = zeros(1, n_heights);
sats_per_plane = zeros(1, n_heights);

fprintf('Computing Walker Star coverage for i = %d deg over %d altitudes...\n', ...
    inclination_deg, n_heights);

for hi = 1:n_heights
    [p, s, ts] = calculate_walker_star_inclined( ...
        heights_km(hi), Min_latitude_deg, Min_elevation_UE, inclination_deg);
    if isinf(ts); ts = NaN; p = NaN; s = NaN; end
    total_sats(hi)     = ts;
    planes(hi)         = p;
    sats_per_plane(hi) = s;
end

%% Save results
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir  = fullfile(fileparts(mfilename('fullpath')), 'simulation_output', ...
                    ['Walker-Star-Inclined_' date_str]);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save(fullfile(out_dir, 'results.mat'), 'total_sats', 'planes', 'sats_per_plane', ...
     'heights_km', 'inclination_deg', 'Min_latitude_deg', 'Min_elevation_UE');
fprintf('Results saved to: %s\n', out_dir);

%% Plot
f1 = figure('Color', 'w', 'Position', [100, 100, 700, 450]);

yyaxis left;
plot(heights_km, planes,         '-', 'Color', '#D95319', 'LineWidth', 1.5); hold on;
plot(heights_km, sats_per_plane, '-', 'Color', '#0072BD', 'LineWidth', 1.5);
ylabel('Planes | Sats/Plane', 'Color', 'k');
ax = gca; ax.YAxis(1).Color = 'k';
ylim([0, max(sats_per_plane, [], 'omitnan') * 1.2]);

yyaxis right;
plot(heights_km, total_sats, '-', 'Color', 'k', 'LineWidth', 2);
ylabel('Total Satellites', 'Color', 'k');
ax.YAxis(2).Color = 'k';
ylim([0, max(total_sats, [], 'omitnan') * 1.5]);

hold off; grid on; ax.GridAlpha = 0.25;
xlabel('Orbital Altitude (km)');
title(sprintf('Walker Star | i = %d\\circ | \\lambda_{min}: %.1f\\circ | \\epsilon_{min}: %.0f\\circ', ...
    inclination_deg, Min_latitude_deg, Min_elevation_UE));
legend('Planes', 'Sats / Plane', 'Total Satellites', 'Location', 'best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

exportgraphics(f1, fullfile(out_dir, 'Walker-Star-Inclined.png'), 'Resolution', 300);
close(f1);
fprintf('Figure saved.\n');
