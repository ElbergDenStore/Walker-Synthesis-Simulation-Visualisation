clear; close all; clc;

% Ensure functions/ subdirectory is on the MATLAB path
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

%% Parameters
heights_km       = 500:0.1:1200;
Min_latitude_deg = 54 + 35/60;   % 54°35'N — southernmost Denmark
Min_elevation_UE = 20;            % degrees minimum elevation angle
inc_array        = 90:-1:84;      % 90°, 89°, 88°, 87°, 86°, 85°, 84°

n_heights = numel(heights_km);
n_inc     = numel(inc_array);

%% Compute results
total_sats_all = zeros(n_heights, n_inc);
fprintf('Computing Walker Star coverage (%d inclinations x %d altitudes)...\n', ...
    n_inc, n_heights);

for ii = 1:n_inc
    for hi = 1:n_heights
        [~, ~, ts] = calculate_walker_star_inclined( ...
            heights_km(hi), Min_latitude_deg, Min_elevation_UE, inc_array(ii));
        if isinf(ts); ts = NaN; end
        total_sats_all(hi, ii) = ts;
    end
    fprintf('  i = %d deg done\n', inc_array(ii));
end

%% Output folder
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir  = fullfile(fileparts(mfilename('fullpath')), 'simulation_output', ...
                    ['Walker-Star-Inclined_' date_str]);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Plotting style
set(0, 'DefaultAxesFontSize', 14);
cmap    = parula(n_inc);
ref_idx = find(inc_array == 90, 1);

%% Figure 1: Total satellites vs altitude
f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 700 450]);
hold on;
h1 = gobjects(n_inc, 1);
for ii = 1:n_inc
    lw       = 1.5 + 1.5 * (inc_array(ii) == 90);   % 90° line is thicker
    h1(ii)   = plot(heights_km, total_sats_all(:, ii), ...
                    'Color', cmap(ii, :), 'LineWidth', lw);
end
hold off;
grid on;
ax1 = gca;
ax1.GridAlpha = 0.25;
xlabel('Orbital Altitude (km)', 'FontName', 'Times New Roman');
ylabel('Satellite Count',       'FontName', 'Times New Roman');
title(sprintf('Walker Star | \\lambda_{min}: %.1f\\circ | \\epsilon_{min}: %.0f\\circ', ...
    Min_latitude_deg, Min_elevation_UE), 'FontName', 'Times New Roman');
labels1 = arrayfun(@(i) sprintf('i = %d\\circ', i), inc_array, 'UniformOutput', false);
lgd1 = legend(h1, labels1, 'Location', 'northeast', ...
              'FontName', 'Times New Roman', 'FontSize', 11);
set(ax1, 'FontName', 'Times New Roman');
exportgraphics(f1, fullfile(out_dir, 'total_satellites.png'), 'Resolution', 300);
close(f1);

%% Figure 2: Satellite penalty vs altitude (relative to i = 90°)
penalty_all = total_sats_all - total_sats_all(:, ref_idx);

f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 700 450]);
hold on;
h2      = gobjects(n_inc - 1, 1);
labels2 = cell(n_inc - 1, 1);
jj = 0;
for ii = 1:n_inc
    if inc_array(ii) == 90; continue; end
    jj          = jj + 1;
    h2(jj)      = plot(heights_km, penalty_all(:, ii), ...
                       'Color', cmap(ii, :), 'LineWidth', 1.5);
    labels2{jj} = sprintf('i = %d\\circ', inc_array(ii));
end
yline(0, 'k--', 'LineWidth', 1.0);
hold off;
grid on;
ax2 = gca;
ax2.GridAlpha = 0.25;
xlabel('Orbital Altitude (km)',          'FontName', 'Times New Roman');
ylabel('Extra Satellites vs i = 90\circ', 'FontName', 'Times New Roman');
title('Walker Star Inclination Penalty',  'FontName', 'Times New Roman');
lgd2 = legend(h2, labels2, 'Location', 'northeast', ...
              'FontName', 'Times New Roman', 'FontSize', 11);
set(ax2, 'FontName', 'Times New Roman');
exportgraphics(f2, fullfile(out_dir, 'satellite_penalty.png'), 'Resolution', 300);
close(f2);

fprintf('Figures saved to: %s\n', out_dir);
