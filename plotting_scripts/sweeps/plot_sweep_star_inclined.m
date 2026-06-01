% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run('plotting_scripts/sweeps/plot_sweep_star_inclined.m')"
close all; clearvars; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'functions'));

%% Output directory
out_dir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'plotting_scripts/figures', 'walker_star_inclined');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Parameters
heights_km       = 500:0.1:1200;
Min_latitude_deg = 54 + 35/60;   % 54°35'N — southernmost Denmark
Min_elevation_UE = 20;            % degrees minimum elevation angle
inc_array        = 90:-1:84;      % 90°, 89°, 88°, 87°, 86°, 85°, 84°

n_heights = numel(heights_km);
n_inc     = numel(inc_array);

%% Compute
total_sats_all = zeros(n_heights, n_inc);
fprintf('Computing Walker Star coverage (%d inclinations x %d altitudes)...\n', n_inc, n_heights);
for ii = 1:n_inc
    for hi = 1:n_heights
        [~, ~, ts] = calculate_walker_star( ...
            heights_km(hi), Min_latitude_deg, Min_elevation_UE, inc_array(ii));
        if isinf(ts); ts = NaN; end
        total_sats_all(hi, ii) = ts;
    end
    fprintf('  i = %d deg done\n', inc_array(ii));
end

%% Plot — multi-inclination overlay
cmap = parula(n_inc);

f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 700 450]);
hold on;
for ii = 1:n_inc
    plot(heights_km, total_sats_all(:, ii), '-', ...
        'Color', cmap(ii, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('i = %d\\circ', inc_array(ii)));
end
hold off;

grid on; set(gca, 'GridAlpha', 0.25);
xlabel('Orbital Altitude (km)', 'FontName', 'Times New Roman');
ylabel('Satellite Count',       'FontName', 'Times New Roman');
title(sprintf('Walker Star | \\lambda_{min}: %.1f\\circ | \\epsilon_{min}: %.0f\\circ', ...
    Min_latitude_deg, Min_elevation_UE), 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 11);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

exportgraphics(f1, fullfile(out_dir, 'walker_star_inclined_sweep.png'), 'Resolution', 300);
close(f1);
fprintf('Saved: walker_star_inclined_sweep.png\n');
fprintf('All figures saved to: %s\n', out_dir);
