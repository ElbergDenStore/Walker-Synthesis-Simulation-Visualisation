% plot_aer_speedup.m
% Benchmarks native MATLAB aer() against custom vectorized ECEF geometry
% across a range of UE counts at a fixed timestep count (2500 steps).
%
% Run headless:
%   xvfb-run -a matlab -nosplash -nodesktop -batch "run('plotting_scripts/plot_aer_speedup.m')"

close all; clearvars; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% =========================================================================
% PARAMETERS
% =========================================================================
height_km   = 1000;
r_earth_m   = 6378.14e3;
min_elev    = 20;       % degrees
sample_time = 60;       % seconds per timestep
num_steps   = 2500;     % 2500 × 60 s ≈ 41.7 h

% Load optimal Walker Delta constellation at this altitude
ws_root = fileparts(fileparts(mfilename('fullpath')));
opt     = load(fullfile(ws_root, 'optimal_constellations.mat'));
[~, ci] = min(abs(opt.heights_km - height_km));
con     = table2struct(opt.best_delta_sats(ci, :));
Total_sats  = con.Sats_per_plane * con.Num_planes;
Num_planes  = con.Num_planes;
Phasing     = con.Phasing;
Inclination = con.Inclination;

start_time = datetime(2025, 6, 1, 12, 0, 0, 'TimeZone', 'UTC');
stop_time  = start_time + seconds((num_steps - 1) * sample_time);

fprintf('=== AER Speedup Benchmark ===\n');
fprintf('Constellation : %d sats | %d planes | %.0f deg inc | %d km\n', ...
    Total_sats, Num_planes, Inclination, height_km);
fprintf('Duration      : %d steps x %d s = %.1f h\n\n', ...
    num_steps, sample_time, num_steps * sample_time / 3600);

% =========================================================================
% UE SWEEP  (num_steps fixed at 2500)
% =========================================================================
ue_counts = [50, 150, 500, 1000, 2500];
n_tests   = numel(ue_counts);

% Pre-generate a pool of UEs covering the Nordic/Arctic region (same as get_cfg)
rng(42);
lat_lims = [54.58, 83.67];   lon_lims = [-60, 30];
n_side   = ceil(sqrt(max(ue_counts) * 2));
[LAT, LON] = meshgrid(linspace(lat_lims(1), lat_lims(2), n_side), ...
                       linspace(lon_lims(1), lon_lims(2), n_side));
pool_lats = LAT(:);  pool_lons = LON(:);
perm      = randperm(numel(pool_lats));   % fixed random draw order

t_aer   = zeros(n_tests, 1);
t_vec   = zeros(n_tests, 1);
t_aer_s = zeros(n_tests, 1);   % aer() : scenario + groundStation setup
t_aer_g = zeros(n_tests, 1);   % aer() : per-UE aer() loop
t_vec_s = zeros(n_tests, 1);   % vec   : scenario + states() prefetch
t_vec_g = zeros(n_tests, 1);   % vec   : per-UE matrix geometry loop
el_err  = zeros(n_tests, 1);
rng_err = zeros(n_tests, 1);

for ti = 1:n_tests
    n_ues   = ue_counts(ti);
    ue_lats = pool_lats(perm(1:n_ues));
    ue_lons = pool_lons(perm(1:n_ues));

    fprintf('[%d/%d] %4d UEs x %d steps ...\n', ti, n_tests, n_ues, num_steps);

    % --- Native aer() ---
    t0 = tic;
      ts = tic;
        sc   = satelliteScenario;
        sc.StartTime  = start_time;
        sc.StopTime   = stop_time;
        sc.SampleTime = sample_time;
        sats  = walkerDelta(sc, height_km*1e3 + r_earth_m, Inclination, ...
                            Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
        ue_gs = groundStation(sc, ue_lats, ue_lons);
      t_aer_s(ti) = toc(ts);
      ts = tic;
        [el_aer, rng_aer] = aer_best(ue_gs, sats, min_elev);
      t_aer_g(ti) = toc(ts);
    t_aer(ti) = toc(t0);
    delete(sc);

    fprintf('  aer()  : %6.2f s  (setup %5.2f s  +  geom %5.2f s)\n', ...
        t_aer(ti), t_aer_s(ti), t_aer_g(ti));

    % --- Vectorized ECEF ---
    t0 = tic;
      ts = tic;
        sc   = satelliteScenario;
        sc.StartTime  = start_time;
        sc.StopTime   = stop_time;
        sc.SampleTime = sample_time;
        sats = walkerDelta(sc, height_km*1e3 + r_earth_m, Inclination, ...
                           Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
        [raw, ~, ~] = states(sats, "CoordinateFrame", "ECEF");
        sat_ecef = permute(raw, [1, 3, 2]);   % 3 x num_sats x nT
      t_vec_s(ti) = toc(ts);
      ts = tic;
        [el_vec, rng_vec] = vec_best(ue_lats, ue_lons, sat_ecef, min_elev);
      t_vec_g(ti) = toc(ts);
    t_vec(ti) = toc(t0);
    delete(sc);

    fprintf('  vec()  : %6.2f s  (setup %5.2f s  +  geom %5.2f s)\n', ...
        t_vec(ti), t_vec_s(ti), t_vec_g(ti));

    % Numerical accuracy vs native aer()
    mask = ~isnan(rng_aer) & ~isnan(rng_vec);
    if any(mask(:))
        rng_err(ti) = max(abs(rng_aer(mask) - rng_vec(mask)));
        el_err(ti)  = max(abs(el_aer(mask)  - el_vec(mask)));
    end

    fprintf('  Speedup: %.1fx  |  max el err: %.2e deg  |  max rng err: %.2e m\n\n', ...
        t_aer(ti)/t_vec(ti), el_err(ti), rng_err(ti));
end

speedup = t_aer ./ t_vec;

% =========================================================================
% SUMMARY TABLE
% =========================================================================
fprintf('%s\n', repmat('=', 1, 86));
fprintf('  %-6s  %8s  %8s  %8s  %8s  %8s  %8s  %9s\n', ...
    'UEs', 'aer (s)', 'aer-setup', 'aer-geom', 'vec (s)', 'vec-setup', 'vec-geom', 'Speedup');
fprintf('%s\n', repmat('-', 1, 86));
for ti = 1:n_tests
    fprintf('  %-6d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.1fx\n', ...
        ue_counts(ti), t_aer(ti), t_aer_s(ti), t_aer_g(ti), ...
        t_vec(ti), t_vec_s(ti), t_vec_g(ti), speedup(ti));
end
fprintf('%s\n', repmat('=', 1, 86));
fprintf('  Max elevation error (deg) :'); fprintf('  %.2e', el_err);  fprintf('\n');
fprintf('  Max range error (m)       :'); fprintf('  %.2e', rng_err); fprintf('\n\n');

% =========================================================================
% FIGURE 1: Computation time vs UE count
% =========================================================================
c_aer = [0.850, 0.325, 0.098];
c_vec = [0,     0.447, 0.741];

f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 580, 380]);
ax1 = axes(f1, 'Color', 'w');

semilogy(ax1, ue_counts, t_aer,   '-o', 'LineWidth', 2.5, 'Color', c_aer, ...
    'MarkerFaceColor', c_aer,   'MarkerSize', 7, 'DisplayName', 'aer()  —  total');
hold(ax1, 'on');
semilogy(ax1, ue_counts, t_aer_g, '--o', 'LineWidth', 1.5, 'Color', c_aer * 0.65, ...
    'MarkerSize', 5, 'DisplayName', 'aer()  —  geometry only');
semilogy(ax1, ue_counts, t_vec,   '-s', 'LineWidth', 2.5, 'Color', c_vec, ...
    'MarkerFaceColor', c_vec,   'MarkerSize', 7, 'DisplayName', 'Vectorized  —  total');
semilogy(ax1, ue_counts, t_vec_g, '--s', 'LineWidth', 1.5, 'Color', c_vec * 0.65, ...
    'MarkerSize', 5, 'DisplayName', 'Vectorized  —  geometry only');

ax1.XScale = 'log';
xlabel(ax1, 'Number of UEs',            'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax1, 'Computation time (s)',      'FontWeight', 'bold', 'FontSize', 13);
legend(ax1, 'Location', 'northwest',    'FontSize', 10);
grid(ax1, 'on');  ax1.Box = 'off';
title(ax1, sprintf('%d sats @ %d km  |  %d timesteps  (\\DeltaT = %d s)', ...
    Total_sats, height_km, num_steps, sample_time), 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f1, fullfile(out_dir, 'aer_speedup_timing.png'), 'Resolution', 300);
close(f1);  fprintf('Saved: aer_speedup_timing.png\n');

% =========================================================================
% FIGURE 2: Speedup factor vs UE count
% =========================================================================
f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 480, 340]);
ax2 = axes(f2, 'Color', 'w');

plot(ax2, ue_counts, speedup, '-^', 'LineWidth', 2.5, ...
    'Color', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880], ...
    'MarkerSize', 8);
yline(ax2, 1, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

ax2.XScale = 'log';
xlabel(ax2, 'Number of UEs',                     'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax2, 'Speedup  (t_{aer} / t_{vec})',       'FontWeight', 'bold', 'FontSize', 13);
grid(ax2, 'on');  ax2.Box = 'off';
title(ax2, sprintf('aer() vs Vectorized  |  %d sats, %d timesteps', Total_sats, num_steps), ...
    'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f2, fullfile(out_dir, 'aer_speedup_factor.png'), 'Resolution', 300);
close(f2);  fprintf('Saved: aer_speedup_factor.png\n');

fprintf('\nAll figures saved to: %s\n', out_dir);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [best_el, best_rng] = aer_best(ue_gs, sats, min_elev)
% Native aer() per UE — returns best-satellite elevation & range (UEs x nT).
    num_ues = numel(ue_gs);
    [~, el0, ~] = aer(ue_gs(1), sats);
    nT = size(el0, 2);
    nS = size(el0, 1);
    best_el  = NaN(num_ues, nT);
    best_rng = NaN(num_ues, nT);
    for i = 1:num_ues
        [~, el_mat, r_mat] = aer(ue_gs(i), sats);
        valid = el_mat >= min_elev;
        r_mat(~valid) = Inf;
        [best_r, bi] = min(r_mat, [], 1);
        has_srv = ~isinf(best_r);
        best_rng(i, has_srv) = best_r(has_srv);
        vc = find(has_srv);
        best_el(i, has_srv) = el_mat(bi(has_srv) + (vc - 1) * nS);
    end
end

function [best_el, best_rng] = vec_best(ue_lats, ue_lons, sat_ecef, min_elev)
% Vectorized ECEF geometry per UE — returns best-satellite elevation & range.
    ue_xyz   = lla2ecef([ue_lats(:), ue_lons(:), zeros(numel(ue_lats), 1)]);
    num_ues  = numel(ue_lats);
    num_sats = size(sat_ecef, 2);
    nT       = size(sat_ecef, 3);
    best_el  = NaN(num_ues, nT);
    best_rng = NaN(num_ues, nT);
    for i = 1:num_ues
        dx   = sat_ecef - ue_xyz(i, :)';          % 3 x num_sats x nT
        slat = sind(ue_lats(i));  clat = cosd(ue_lats(i));
        slon = sind(ue_lons(i));  clon = cosd(ue_lons(i));
        R_enu = [-slon,           clon,          0; ...
                 -slat*clon,     -slat*slon,     clat; ...
                  clat*clon,      clat*slon,     slat];
        ve   = reshape(R_enu * reshape(dx, 3, []), 3, num_sats, nT);
        E    = reshape(ve(1,:,:), num_sats, nT);
        N    = reshape(ve(2,:,:), num_sats, nT);
        U    = reshape(ve(3,:,:), num_sats, nT);
        r_mat  = sqrt(E.^2 + N.^2 + U.^2);
        el_mat = asind(U ./ r_mat);
        valid  = el_mat >= min_elev;
        r_mat(~valid) = Inf;
        [best_r, bi] = min(r_mat, [], 1);
        has_srv = ~isinf(best_r);
        best_rng(i, has_srv) = best_r(has_srv);
        vc = find(has_srv);
        best_el(i, has_srv) = el_mat(bi(has_srv) + (vc - 1) * num_sats);
    end
end
