% TEMP_COMPARE_ANALYTICAL_FORMULAS
% Temporary script: overlay three series on one plot:
%   (1) Analytical from calculate_walker_star          (saved in sweep .mat)
%   (2) Analytical from calculate_walker_star_inclined at i=90 (recomputed)
%   (3) Numerical gridsearch best
%
% Delete this file when done.

addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

workspace_root = fileparts(mfilename('fullpath'));
sim_dir        = fullfile(workspace_root, 'simulation_output');

%% ---- Locate most recent Walker Star sweep ------------------------------
hits = dir(fullfile(sim_dir, 'Master_Sweep_WalkerStar_*'));
hits = hits([hits.isdir]);
if isempty(hits)
    error('No Master_Sweep_WalkerStar_* folders found in:\n  %s', sim_dir);
end
[~, idx]     = max([hits.datenum]);
sweep_folder = fullfile(sim_dir, hits(idx).name);
fprintf('Sweep folder: %s\n', sweep_folder);

loaded = load(fullfile(sweep_folder, 'Master_Altitude_Sweep_Results.mat'));

heights_km     = loaded.heights_km;
star_sats      = loaded.star_sats;
best_star_sats = loaded.best_star_sats;

[heights_km, sort_idx] = sort(heights_km, 'ascend');
star_sats      = star_sats(sort_idx);
best_star_sats = best_star_sats(sort_idx, :);

%% ---- Extract config params for recomputing inclined formula ------------
cfg      = loaded.Master_config;
min_lat  = min(cfg.Lat_range_deg);
min_elev = cfg.Min_elevation_UE;

%% ---- Recompute inclined@90 for every altitude --------------------------
n          = numel(heights_km);
incl90_T   = zeros(n, 1);
for k = 1:n
    [~, ~, incl90_T(k)] = calculate_walker_star( ...
        heights_km(k), min_lat, min_elev, 90);
end

anal_T = [star_sats.Total_sats]';
num_T  = best_star_sats.Total_sats;

%% ---- Plot --------------------------------------------------------------
title_str = sprintf( ...
    'Walker Star Analytical Comparison  |  \\lambda: %.1f°  |  \\epsilon_{min}: %.1f°', ...
    min_lat, min_elev);

set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

figure('Color', 'w', 'Position', [100 100 760 480]);
hold on;
plot(heights_km, anal_T,   'r-o', 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Analytical (calculate\_walker\_star)');
plot(heights_km, incl90_T, 'm-s', 'MarkerFaceColor', 'm', ...
    'DisplayName', 'Analytical (calculate\_walker\_star\_inclined, i=90°)');
plot(heights_km, num_T,    'b-^', 'MarkerFaceColor', 'b', ...
    'DisplayName', 'Numerical gridsearch');
xlabel('Equatorial Orbital Altitude (km)');
ylabel('Satellite Count');
title(title_str);
legend('Location', 'northeast');
grid on; hold off;

plot_path = fullfile(sweep_folder, 'Star_SMAD_vs_Inclined90_vs_Numerical.png');
exportgraphics(gcf, plot_path, 'Resolution', 300);
fprintf('Plot saved to: %s\n', plot_path);

%% ---- Console diff table ------------------------------------------------
fprintf('\n%-8s  %-10s  %-10s  %-10s  %-12s  %-12s\n', ...
    'Alt(km)', 'SMAD(T)', 'Incl90(T)', 'Num(T)', ...
    'SMAD-Num', 'Incl90-Num');
fprintf('%s\n', repmat('-', 1, 68));
for k = 1:n
    fprintf('%-8d  %-10d  %-10d  %-10d  %-12d  %-12d\n', ...
        heights_km(k), anal_T(k), incl90_T(k), num_T(k), ...
        anal_T(k) - num_T(k), incl90_T(k) - num_T(k));
end
