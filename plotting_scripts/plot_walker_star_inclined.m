% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run('plotting_scripts/plot_walker_star_inclined.m')"
close all; clearvars; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

%% Output Directory
out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures', 'walker_star_inclined');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% Parameters
heights_km       = 500:0.1:1200;
Min_latitude_deg = 54 + 35/60;   % 54°35'N — southernmost Denmark
Min_elevation_UE = 20;
inclination_deg  = 90;            % single inclination to evaluate

n_heights = numel(heights_km);

%% Compute
total_sats   = zeros(1, n_heights);
planes       = zeros(1, n_heights);
sats_per_plane = zeros(1, n_heights);

fprintf('Computing Walker Star coverage for i = %d deg over %d altitudes...\n', inclination_deg, n_heights);
for hi = 1:n_heights
    [p, s, ts] = calculate_walker_star_inclined( ...
        heights_km(hi), Min_latitude_deg, Min_elevation_UE, inclination_deg);
    if isinf(ts); ts = NaN; p = NaN; s = NaN; end
    total_sats(hi)    = ts;
    planes(hi)        = p;
    sats_per_plane(hi) = s;
end

%% Plotting style
set(0, 'DefaultAxesFontSize', 14);

%% Figure 1: Planes, Sats/Plane, and Total Satellites vs altitude
f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 700 450]);

yyaxis left;
plot(heights_km, planes,         '-', 'Color', '#D95319', 'LineWidth', 1.5); hold on;
plot(heights_km, sats_per_plane, '-', 'Color', '#0072BD', 'LineWidth', 1.5);
ylabel('Planes | Sats/Plane', 'Color', 'k', 'FontName', 'Times New Roman');
ax = gca;
ax.YAxis(1).Color = 'k';
ylim([0, max(sats_per_plane, [], 'omitnan') * 1.2]);

yyaxis right;
plot(heights_km, total_sats, '-', 'Color', 'k', 'LineWidth', 2);
ylabel('Total Satellites', 'Color', 'k', 'FontName', 'Times New Roman');
ax.YAxis(2).Color = 'k';
ylim([0, max(total_sats, [], 'omitnan') * 1.5]);

hold off; grid on; ax.GridAlpha = 0.25;
xlabel('Orbital Altitude (km)', 'FontName', 'Times New Roman');
title(sprintf('Walker Star | i = %d\\circ | \\lambda_{min}: %.1f\\circ | \\epsilon_{min}: %.0f\\circ', ...
    inclination_deg, Min_latitude_deg, Min_elevation_UE), 'FontName', 'Times New Roman');
legend('Planes', 'Sats / Plane', 'Total Satellites', 'Location', 'best', ...
       'FontName', 'Times New Roman', 'FontSize', 11);
set(gca, 'FontName', 'Times New Roman');
exportgraphics(f1, fullfile(out_dir, 'walker_star_TPS.png'), 'Resolution', 300);
close(f1);
fprintf('Saved: walker_star_TPS.png\n');

%% Figure 2: Total satellites vs altitude (clean single-line view)
f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 500 300]);
plot(heights_km, total_sats, '-', 'Color', '#0072BD', 'LineWidth', 2);
grid on; gca().GridAlpha = 0.25;
xlabel('Orbital Altitude (km)', 'FontName', 'Times New Roman');
ylabel('Total Satellites',      'FontName', 'Times New Roman');
title(sprintf('Walker Star | i = %d\\circ | \\lambda_{min}: %.1f\\circ | \\epsilon_{min}: %.0f\\circ', ...
    inclination_deg, Min_latitude_deg, Min_elevation_UE), 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');
exportgraphics(f2, fullfile(out_dir, 'total_satellites.png'), 'Resolution', 300);
close(f2);
fprintf('Saved: total_satellites.png\n');

fprintf('All figures saved to: %s\n', out_dir);
