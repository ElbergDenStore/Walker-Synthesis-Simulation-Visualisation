% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run('plotting_scripts/plot_tx_gain_vs_altitude.m')"
close all; clearvars; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

%% Output directory
out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'plotting_scripts/figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Parameters
heights_km = 500:1:1200;
elev_nadir = 90;    % nadir: slant range = altitude
n_h        = numel(heights_km);

% Use reference frequency (20 GHz) — frequency only shifts gain vertically by a constant
freq_Hz    = 20e9;
G_elem_dBi = 6;    % assumed element gain for array size estimate

%% Compute gain
G = zeros(1, n_h);
for hi = 1:n_h
    G(hi) = get_adjusted_tx_gain(heights_km(hi) * 1e3, elev_nadir, freq_Hz);
end

%% Primary axes — Tx Gain
f1  = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 700 450]);
ax1 = axes(f1, 'Position', [0.11 0.12 0.67 0.74]);   % leave right margin for ax2 label
plot(ax1, heights_km, G, '-', 'Color', '#0072BD', 'LineWidth', 2);
grid(ax1, 'on');
set(ax1, 'GridAlpha', 0.25, 'FontName', 'Times New Roman', 'FontSize', 14, 'Box', 'off');
xlabel(ax1, 'Orbital Altitude (km)',   'FontName', 'Times New Roman');
ylabel(ax1, 'Required Tx Gain (dBi)', 'FontName', 'Times New Roman');
title(ax1,  'Tx Gain to Maintain Constant Nadir Beam Footprint', 'FontName', 'Times New Roman');

%% Derived right axis — element count (log scale)
% N = 10^((G - G_elem)/10)  →  log10(N) is linear in G
% Setting log-scale limits to the transformed left-axis limits ensures perfect alignment.
yl    = ylim(ax1);
N_lim = 10.^((yl - G_elem_dBi) / 10);

ax2 = axes(f1, 'Position', ax1.Position, ...
               'YAxisLocation', 'right', ...
               'Color', 'none', 'YScale', 'log', ...
               'XTick', [], 'YLim', N_lim, 'XLim', ax1.XLim, ...
               'Box', 'off');
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel(ax2, 'Approximate Array Elements', 'FontName', 'Times New Roman');

% Doubling ticks starting at 500: 500, 1000, 2000, 4000, ...
base_ticks = [1000, 2000, 3000]
base_ticks = base_ticks(base_ticks >= N_lim(1) & base_ticks <= N_lim(2));
if ~isempty(base_ticks)
    set(ax2, 'YTick', base_ticks, ...
             'YTickLabel', arrayfun(@num2str, base_ticks, 'UniformOutput', false));
end

%% Export
exportgraphics(f1, fullfile(out_dir, 'tx_gain_vs_altitude.png'), 'Resolution', 300);
close(f1);
fprintf('Saved: tx_gain_vs_altitude.png\n');
