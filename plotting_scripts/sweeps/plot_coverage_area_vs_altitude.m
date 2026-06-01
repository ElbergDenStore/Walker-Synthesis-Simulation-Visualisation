% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run('plotting_scripts/sweeps/plot_coverage_area_vs_altitude.m')"
close all; clearvars; clc;

%% Output directory
out_dir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'plotting_scripts/figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Parameters
heights_km  = 500:1:1200;
elev_min    = 20;          % minimum elevation angle (degrees)
Re_km       = 6371;        % Earth radius (km)

%% Compute spherical cap coverage area
% Step 1: nadir angle η at the satellite
eta_rad  = asin(Re_km * cosd(elev_min) ./ (Re_km + heights_km));

% Step 2: Earth central half-angle ρ = 90° - ε - η
rho_rad  = deg2rad(90 - elev_min) - eta_rad;

% Spherical cap area (km²)
A_km2    = 2 * pi * Re_km^2 .* (1 - cos(rho_rad));

%% Primary axes — coverage area (LEO sweep)
f1  = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 700 450]);
ax1 = axes(f1, 'Position', [0.11 0.12 0.78 0.74]);
plot(ax1, heights_km, A_km2 / 1e6, '-', 'Color', '#0072BD', 'LineWidth', 2);
grid(ax1, 'on');
set(ax1, 'GridAlpha', 0.25, 'FontName', 'Times New Roman', 'FontSize', 14, 'Box', 'off');
xlabel(ax1, 'Orbital Altitude (km)',     'FontName', 'Times New Roman');
ylabel(ax1, 'Coverage Area (10^6 km^2)', 'FontName', 'Times New Roman');
title(ax1,  sprintf('Single-Satellite Coverage | \\epsilon_{min} = %d\\circ', elev_min), ...
      'FontName', 'Times New Roman');

exportgraphics(f1, fullfile(out_dir, 'coverage_area_vs_altitude.png'), 'Resolution', 300);
close(f1);
fprintf('Saved: coverage_area_vs_altitude.png\n');

%% Figure 2 — full orbital regime (1 km to 50 000 km)
h_wide_km   = 1:1:50000;
eta_wide    = asin(Re_km * cosd(elev_min) ./ (Re_km + h_wide_km));
rho_wide    = deg2rad(90 - elev_min) - eta_wide;
A_wide_km2  = 2 * pi * Re_km^2 .* (1 - cos(rho_wide));

f2  = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 700 450]);
ax3 = axes(f2, 'Position', [0.11 0.12 0.78 0.74]);
plot(ax3, h_wide_km, A_wide_km2 / 1e6, '-', 'Color', '#0072BD', 'LineWidth', 2);
grid(ax3, 'on');
set(ax3, 'GridAlpha', 0.25, 'FontName', 'Times New Roman', 'FontSize', 14, 'Box', 'off');
xlabel(ax3, 'Orbital Altitude (km)',     'FontName', 'Times New Roman');
ylabel(ax3, 'Coverage Area (10^6 km^2)', 'FontName', 'Times New Roman');
title(ax3,  sprintf('Single-Satellite Coverage | \\epsilon_{min} = %d\\circ', elev_min), ...
      'FontName', 'Times New Roman');

% Regime tick marks and vertical lines
regimes     = [1000,  20000,  35786];
regime_lbls = {'LEO', 'MEO',  'GEO'};
colors_reg  = {'#A2142F', '#EDB120', '#77AC30'};
yl3 = ylim(ax3);
hold(ax3, 'on');
for ri = 1:numel(regimes)
    xline(ax3, regimes(ri), '--', 'Color', colors_reg{ri}, 'LineWidth', 1.2, ...
          'Label', regime_lbls{ri}, 'LabelVerticalAlignment', 'bottom', ...
          'FontName', 'Times New Roman', 'FontSize', 12);
end
hold(ax3, 'off');

exportgraphics(f2, fullfile(out_dir, 'coverage_area_vs_altitude_wide.png'), 'Resolution', 300);
close(f2);
fprintf('Saved: coverage_area_vs_altitude_wide.png\n');
