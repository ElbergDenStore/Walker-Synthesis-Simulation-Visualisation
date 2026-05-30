% plot_aer_speedup.m
% Benchmarks four AER computation methods across UE counts at fixed 2500 timesteps,
% plus a parallel gridsearch case study: 32 independent 50-UE simulations.
%
% Methods:
%   1  aer()         — native MATLAB toolbox aer(), per-UE loop
%   2  vec-toolbox   — toolbox states() ECEF prefetch + vectorised ENU loop
%   3  vec-fastmath  — fast_walker_ecef() (zero toolbox overhead) + same loop
%   4  par-gridsearch — 32 × 50-UE fastmath sims dispatched via parfor
%                       (models constellation gridsearch: one config per worker)
%
% Run headless:
%   xvfb-run -a matlab -nosplash -nodesktop -batch "run('plotting_scripts/plot_aer_speedup.m')"

close all; clearvars; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'plotting_scripts/figures/optimisations');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% =========================================================================
% PARAMETERS
% =========================================================================
height_km   = 1000;
r_earth_m   = 6378.14e3;
min_elev    = 20;
sample_time = 60;       % seconds per timestep
num_steps   = 2500;

ue_counts = [50, 150, 500, 1000, 2500];
n_tests   = numel(ue_counts);

% =========================================================================
% LOAD CONSTELLATION
% =========================================================================
ws_root = fileparts(fileparts(mfilename('fullpath')));
opt     = load(fullfile(ws_root, 'optimal_constellations.mat'));
[~, ci] = min(abs(opt.heights_km - height_km));
con     = table2struct(opt.best_delta_sats(ci, :));

Total_sats     = con.Sats_per_plane * con.Num_planes;
Num_planes     = con.Num_planes;
Sats_per_plane = con.Sats_per_plane;
Phasing        = con.Phasing;
Inclination    = con.Inclination;

start_time     = datetime(2025, 6, 1, 12, 0, 0, 'TimeZone', 'UTC');
stop_time      = start_time + seconds((num_steps - 1) * sample_time);
time_steps_sec = 0 : sample_time : (num_steps - 1) * sample_time;

fprintf('=== AER Method Benchmark ===\n');
fprintf('Constellation: %d sats | %d planes | %.0f deg inc | %d km\n', ...
    Total_sats, Num_planes, Inclination, height_km);
fprintf('Timesteps    : %d x %d s = %.1f h\n\n', ...
    num_steps, sample_time, num_steps * sample_time / 3600);

% =========================================================================
% UE POOL  (fixed draw order for reproducibility)
% =========================================================================
rng(42);
lat_lims = [54.58, 83.67];  lon_lims = [-60, 30];
n_side   = ceil(sqrt(max(ue_counts) * 2));
[LAT, LON] = meshgrid(linspace(lat_lims(1), lat_lims(2), n_side), ...
                       linspace(lon_lims(1), lon_lims(2), n_side));
pool_lats = LAT(:);  pool_lons = LON(:);
perm      = randperm(numel(pool_lats));

% =========================================================================
% PRE-COMPUTE SATELLITE POSITIONS via fast_walker_ecef  (methods 3 & 4)
% Run once — independent of UE count, so a single timing is fair.
% =========================================================================
fprintf('Pre-computing satellite positions via fast_walker_ecef...\n');
fast_setup_tic = tic;
sat_ecef_fast = fast_walker_ecef(height_km*1e3, Inclination, Num_planes, ...
                                  Sats_per_plane, Phasing, time_steps_sec, start_time);
t_fast_setup = toc(fast_setup_tic);
fprintf('  Done: %.2f s  (%d sats x %d steps)\n\n', ...
    t_fast_setup, size(sat_ecef_fast,2), size(sat_ecef_fast,3));

% =========================================================================
% SWEEP OVER UE COUNTS  (methods 1-3)
% =========================================================================
% t_setup(m, ti) = orbit/scenario setup time for method m at test ti
% t_geom(m, ti)  = geometry loop time
% t_total(m, ti) = end-to-end time (setup + geom)
t_setup = zeros(3, n_tests);
t_geom  = zeros(3, n_tests);
t_total = zeros(3, n_tests);

for ti = 1:n_tests
    n_ues   = ue_counts(ti);
    ue_lats = pool_lats(perm(1:n_ues));
    ue_lons = pool_lons(perm(1:n_ues));

    fprintf('[%d/%d]  %4d UEs x %d steps\n', ti, n_tests, n_ues, num_steps);

    % ------------------------------------------------------------------
    % Method 1: native aer()
    %   Setup: satelliteScenario + walkerDelta + groundStation
    %   Geom:  aer() called once per UE
    % ------------------------------------------------------------------
    t0 = tic;
      ts = tic;
        sc1 = satelliteScenario;
        sc1.StartTime = start_time;  sc1.StopTime = stop_time;  sc1.SampleTime = sample_time;
        sats1  = walkerDelta(sc1, height_km*1e3 + r_earth_m, Inclination, ...
                             Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
        ue_gs  = groundStation(sc1, ue_lats, ue_lons);
      t_setup(1, ti) = toc(ts);
      ts = tic;
        run_aer(ue_gs, sats1, min_elev);
      t_geom(1, ti) = toc(ts);
    t_total(1, ti) = toc(t0);
    delete(sc1);

    % ------------------------------------------------------------------
    % Method 2: vec-toolbox
    %   Setup: satelliteScenario + walkerDelta + states() ECEF prefetch
    %   Geom:  vectorised ENU loop (no toolbox calls)
    % ------------------------------------------------------------------
    t0 = tic;
      ts = tic;
        sc2 = satelliteScenario;
        sc2.StartTime = start_time;  sc2.StopTime = stop_time;  sc2.SampleTime = sample_time;
        sats2 = walkerDelta(sc2, height_km*1e3 + r_earth_m, Inclination, ...
                            Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
        [raw, ~, ~]   = states(sats2, "CoordinateFrame", "ECEF");
        sat_ecef_tb   = permute(raw, [1, 3, 2]);
      t_setup(2, ti) = toc(ts);
      ts = tic;
        run_vec(ue_lats, ue_lons, sat_ecef_tb, min_elev);
      t_geom(2, ti) = toc(ts);
    t_total(2, ti) = toc(t0);
    delete(sc2);

    % ------------------------------------------------------------------
    % Method 3: vec-fastmath
    %   Setup: fast_walker_ecef()  (pre-computed once above, reused here)
    %   Geom:  same vectorised ENU loop as method 2
    % ------------------------------------------------------------------
    t_setup(3, ti) = t_fast_setup;   % orbit generation, independent of UE count
    ts = tic;
      run_vec(ue_lats, ue_lons, sat_ecef_fast, min_elev);
    t_geom(3, ti)  = toc(ts);
    t_total(3, ti) = t_fast_setup + t_geom(3, ti);

    fprintf('  aer():        %6.2f s  (setup %5.2f  +  geom %5.2f)\n', t_total(1,ti), t_setup(1,ti), t_geom(1,ti));
    fprintf('  vec-toolbox:  %6.2f s  (setup %5.2f  +  geom %5.2f)\n', t_total(2,ti), t_setup(2,ti), t_geom(2,ti));
    fprintf('  vec-fastmath: %6.2f s  (setup %5.2f  +  geom %5.2f)\n\n', t_total(3,ti), t_fast_setup, t_geom(3,ti));
end

% =========================================================================
% METHOD 4: par-gridsearch
%   Models the gridsearch use case: many independent constellation configs
%   evaluated in parallel, each a full fast_walker_ecef + ENU geometry sim.
%
%   n_per_worker sims are assigned per worker so that per-iteration work
%   (n_per_worker * ~0.1 s) dominates worker broadcast/dispatch overhead.
%   sat_ecef_fast is broadcast once to all workers at parfor launch.
%
%   NOTE: with only 32 sims (1 per worker, ~0.12 s each) overhead dominates
%   and speedup ≈ 1×.  Running n_per_worker=10 sims per worker gives enough
%   work (~1.2 s/worker) for the parallel benefit to show clearly.
% =========================================================================
n_ues_cs     = 50;    % UEs per individual simulation
n_per_worker = 10;    % simulations assigned to each worker
n_sims       = n_per_worker * 32;   % total = 320

% Build per-simulation UE subsets (different slice of the pool each time)
ue_lats_cs = cell(n_sims, 1);
ue_lons_cs = cell(n_sims, 1);
for k = 1:n_sims
    i0 = mod((k-1) * n_ues_cs, numel(pool_lats) - n_ues_cs) + 1;
    ue_lats_cs{k} = pool_lats(perm(i0 : i0+n_ues_cs-1));
    ue_lons_cs{k} = pool_lons(perm(i0 : i0+n_ues_cs-1));
end

% --- Sequential baseline: all sims one after another ---
fprintf('Method 4 — sequential baseline: %d x %d-UE fastmath sims...\n', n_sims, n_ues_cs);
t_seq_tic = tic;
for k = 1:n_sims
    run_vec(ue_lats_cs{k}, ue_lons_cs{k}, sat_ecef_fast, min_elev);
end
t_seq_total = toc(t_seq_tic);
fprintf('  sequential: %.2f s  (%.3f s / sim)\n\n', t_seq_total, t_seq_total/n_sims);

% --- Parallel: all sims dispatched over workers ---
fprintf('Starting parallel pool...\n');
if isempty(gcp('nocreate'))
    parpool();
end
% Warm-up: two iterations so scheduler spin-up is not counted in timing
parfor k = 1:2
    run_vec(ue_lats_cs{k}, ue_lons_cs{k}, sat_ecef_fast, min_elev);
end

fprintf('Method 4 — parallel gridsearch: %d x %d-UE fastmath sims (parfor, ~%d/worker)...\n', ...
    n_sims, n_ues_cs, n_per_worker);
t_par_tic = tic;
parfor k = 1:n_sims
    run_vec(ue_lats_cs{k}, ue_lons_cs{k}, sat_ecef_fast, min_elev);
end
t_par_total_cs = toc(t_par_tic);
fprintf('  parallel:   %.2f s  (%.3f s / sim)\n', t_par_total_cs, t_par_total_cs/n_sims);
par_speedup = t_seq_total / t_par_total_cs;
fprintf('  Parallel speedup: %.1fx over sequential\n\n', par_speedup);

% =========================================================================
% METHOD 5: par-toolbox  (2500 UEs x 2500 steps)
%   Direct equivalent of coverage_simulator_function.m with use_parallel=true.
%   Setup: satelliteScenario + walkerDelta + states() ECEF prefetch  (same as
%          method 2, not accelerated by parallelism).
%   Geom:  parfor over UEs — identical ENU math to run_vec, dispatched over
%          the same worker pool (parfor(..., Inf) lets MATLAB auto-assign).
% =========================================================================
n_ues_par5  = 2500;
ue_lats_par5 = pool_lats(perm(1:n_ues_par5));
ue_lons_par5 = pool_lons(perm(1:n_ues_par5));

fprintf('Method 5 — par-toolbox, %d UEs x %d steps...\n', n_ues_par5, num_steps);
t0_par_tb = tic;
  ts = tic;
    sc5   = satelliteScenario;
    sc5.StartTime = start_time;  sc5.StopTime = stop_time;  sc5.SampleTime = sample_time;
    sats5 = walkerDelta(sc5, height_km*1e3 + r_earth_m, Inclination, ...
                        Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
    [raw5, ~, ~]    = states(sats5, "CoordinateFrame", "ECEF");
    sat_ecef_par_tb = permute(raw5, [1, 3, 2]);
  t_par_tb_setup = toc(ts);
  ts = tic;
    run_vec_par_ues(ue_lats_par5, ue_lons_par5, sat_ecef_par_tb, min_elev);
  t_par_tb_geom  = toc(ts);
t_par_tb_total = toc(t0_par_tb);
delete(sc5);

tb_seq_total = t_setup(2,end) + t_geom(2,end);   % vec-toolbox at 2500 UEs
tb_speedup   = tb_seq_total / t_par_tb_total;
fprintf('  par-toolbox: %.2f s  (setup %5.2f  +  geom %5.2f)\n', t_par_tb_total, t_par_tb_setup, t_par_tb_geom);
fprintf('  vs vec-toolbox: %.2f s  ->  Speedup: %.1fx\n\n', tb_seq_total, tb_speedup);

% =========================================================================
% SUMMARY TABLE
% =========================================================================
method_names_sweep = {'aer()', 'vec-toolbox', 'vec-fastmath'};
W = 78;
fprintf('%s\n', repmat('=', 1, W));
fprintf('  %-14s  %6s  %10s  %10s  %10s  %9s\n', ...
    'Method', 'UEs', 'Setup (s)', 'Geom (s)', 'Total (s)', 'Speedup');
fprintf('%s\n', repmat('-', 1, W));
for ti = 1:n_tests
    for m = 1:3
        sp = t_total(1,ti) / t_total(m,ti);
        fprintf('  %-14s  %6d  %10.2f  %10.2f  %10.2f  %8.1fx\n', ...
            method_names_sweep{m}, ue_counts(ti), t_setup(m,ti), t_geom(m,ti), t_total(m,ti), sp);
    end
    fprintf('%s\n', repmat('-', 1, W));
end
fprintf('%s\n', repmat('=', 1, W));
fprintf('\n  par-toolbox case study: %d UEs x %d steps\n', n_ues_par5, num_steps);
fprintf('  vec-toolbox (seq): %.2f s  |  par-toolbox: %.2f s  |  Speedup: %.1fx\n', ...
    tb_seq_total, t_par_tb_total, tb_speedup);
fprintf('\n  Gridsearch case study: %d independent %d-UE fastmath sims\n', n_sims, n_ues_cs);
fprintf('  Sequential : %.2f s  |  Parallel: %.2f s  |  Speedup: %.1fx\n', ...
    t_seq_total, t_par_total_cs, par_speedup);
fprintf('%s\n\n', repmat('=', 1, W));

speedup_tb   = t_total(1,:) ./ t_total(2,:);
speedup_fast = t_total(1,:) ./ t_total(3,:);

% =========================================================================
% FIGURE 1: Computation time vs UE count (log-log)
% =========================================================================
c1 = [0.850, 0.325, 0.098];
c2 = [0,     0.447, 0.741];
c3 = [0.4660, 0.6740, 0.1880];

f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 620, 420]);
ax1 = axes(f1, 'Color', 'w');
hold(ax1, 'on');

semilogy(ax1, ue_counts, t_total(1,:), '-o',  'LineWidth', 2.5, 'Color', c1,      'MarkerFaceColor', c1,      'MarkerSize', 7, 'DisplayName', 'aer()  — total');
semilogy(ax1, ue_counts, t_geom(1,:),  '--o', 'LineWidth', 1.5, 'Color', c1*0.55, 'MarkerSize', 5,            'DisplayName', 'aer()  — geometry only');
semilogy(ax1, ue_counts, t_total(2,:), '-s',  'LineWidth', 2.5, 'Color', c2,      'MarkerFaceColor', c2,      'MarkerSize', 7, 'DisplayName', 'vec-toolbox  — total');
semilogy(ax1, ue_counts, t_geom(2,:),  '--s', 'LineWidth', 1.5, 'Color', c2*0.55, 'MarkerSize', 5,            'DisplayName', 'vec-toolbox  — geometry only');
semilogy(ax1, ue_counts, t_total(3,:), '-^',  'LineWidth', 2.5, 'Color', c3,      'MarkerFaceColor', c3,      'MarkerSize', 7, 'DisplayName', 'vec-fastmath  — total');
semilogy(ax1, ue_counts, t_geom(3,:),  '--^', 'LineWidth', 1.5, 'Color', c3*0.55, 'MarkerSize', 5,            'DisplayName', 'vec-fastmath  — geometry only');

ax1.XScale = 'log';
xlabel(ax1, 'Number of UEs',         'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax1, 'Computation time (s)',   'FontWeight', 'bold', 'FontSize', 13);
legend(ax1, 'Location', 'northwest', 'FontSize', 9);
grid(ax1, 'on');  ax1.Box = 'off';
title(ax1, sprintf('%d sats @ %d km  |  %d timesteps  (\\DeltaT = %d s)', ...
    Total_sats, height_km, num_steps, sample_time), 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f1, fullfile(out_dir, 'aer_speedup_timing.png'), 'Resolution', 300);
close(f1);  fprintf('Saved: aer_speedup_timing.png\n');

% =========================================================================
% FIGURE 2: Speedup factor vs UE count
% =========================================================================
f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 500, 360]);
ax2 = axes(f2, 'Color', 'w');
hold(ax2, 'on');

plot(ax2, ue_counts, speedup_tb,   '-s', 'LineWidth', 2.5, 'Color', c2, 'MarkerFaceColor', c2, 'MarkerSize', 8, 'DisplayName', 'vec-toolbox');
plot(ax2, ue_counts, speedup_fast, '-^', 'LineWidth', 2.5, 'Color', c3, 'MarkerFaceColor', c3, 'MarkerSize', 8, 'DisplayName', 'vec-fastmath');
yline(ax2, 1, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

ax2.XScale = 'log';
xlabel(ax2, 'Number of UEs',                    'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax2, 'Speedup  (t_{aer} / t_{method})',   'FontWeight', 'bold', 'FontSize', 13);
legend(ax2, 'Location', 'northwest', 'FontSize', 11);
grid(ax2, 'on');  ax2.Box = 'off';
title(ax2, sprintf('Speedup relative to aer()  |  %d sats, %d timesteps', Total_sats, num_steps), ...
    'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f2, fullfile(out_dir, 'aer_speedup_factor.png'), 'Resolution', 300);
close(f2);  fprintf('Saved: aer_speedup_factor.png\n');

% =========================================================================
% FIGURE 3: Case study — two types of parallelism (4 bars)
%   Left pair:  single sim 2500×2500  —  vec-toolbox vs par-toolbox
%   Right pair: gridsearch 32×50      —  seq-fastmath vs par-fastmath
% =========================================================================
% Bar data: [setup_time, geom_time] for each of 4 scenarios
bar_setup4 = [t_setup(2,end),   t_par_tb_setup,   t_fast_setup,  t_fast_setup];
bar_geom4  = [t_geom(2,end),    t_par_tb_geom,    t_seq_total,   t_par_total_cs];

bar_labels = {'vec-toolbox', 'par-toolbox', ...
              sprintf('seq (%d\\times%d UE)', n_sims, n_ues_cs), ...
              sprintf('par (%d\\times%d UE)', n_sims, n_ues_cs)};

f3 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 600, 390]);
ax3 = axes(f3, 'Color', 'w');

bh = bar(ax3, [bar_setup4; bar_geom4]', 'stacked');
bh(1).FaceColor = [0.72, 0.72, 0.72];  bh(1).DisplayName = 'Setup / orbit generation';
bh(2).FaceColor = [0.18, 0.38, 0.60];  bh(2).DisplayName = 'Geometry (total)';

% Annotate speedups above each bar
totals4  = bar_setup4 + bar_geom4;
max_bar4 = max(totals4);
text(ax3, 1, totals4(1) + max_bar4*0.03, '1\times',               'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
text(ax3, 2, totals4(2) + max_bar4*0.03, sprintf('%.0f\times', tb_speedup),  'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', c2);
text(ax3, 3, totals4(3) + max_bar4*0.03, '1\times',               'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
text(ax3, 4, totals4(4) + max_bar4*0.03, sprintf('%.0f\times', par_speedup), 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', c3);

% Group separator line
xl = xline(ax3, 2.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');

ax3.XTickLabel         = bar_labels;
ax3.XTickLabelRotation = 8;
ylabel(ax3, 'Total wall-clock time (s)', 'FontWeight', 'bold', 'FontSize', 13);
legend(ax3, 'Location', 'northeast', 'FontSize', 11);
grid(ax3, 'on');  ax3.Box = 'off';  ax3.XGrid = 'off';
title(ax3, sprintf('%d sats @ %d km  |  left: 2500 UEs \\times %d steps  |  right: %d\\times%d-UE sims (~%d/worker)', ...
    Total_sats, height_km, num_steps, n_sims, n_ues_cs, n_per_worker), 'FontSize', 10, 'FontWeight', 'bold');

exportgraphics(f3, fullfile(out_dir, 'aer_speedup_case_study.png'), 'Resolution', 300);
close(f3);  fprintf('Saved: aer_speedup_case_study.png\n');

fprintf('\nAll figures saved to: %s\n', out_dir);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function run_vec_par_ues(ue_lats, ue_lons, sat_ecef, min_elev)
% Vectorised ECEF geometry, parfor over UEs — mirrors coverage_simulator_function.m
% with use_parallel=true.  Identical ENU math to run_vec; parallelism across UEs.
    ue_xyz   = lla2ecef([ue_lats(:), ue_lons(:), zeros(numel(ue_lats), 1)]);
    num_ues  = numel(ue_lats);
    num_sats = size(sat_ecef, 2);
    nT       = size(sat_ecef, 3);
    best_el  = NaN(num_ues, nT);
    best_rng = NaN(num_ues, nT);
    parfor (i = 1:num_ues, Inf)
        dx    = sat_ecef - ue_xyz(i, :)';
        slat  = sind(ue_lats(i));  clat = cosd(ue_lats(i));
        slon  = sind(ue_lons(i));  clon = cosd(ue_lons(i));
        R_enu = [-slon,       clon,       0; ...
                 -slat*clon, -slat*slon,  clat; ...
                  clat*clon,  clat*slon,  slat];
        ve    = reshape(R_enu * reshape(dx, 3, []), 3, num_sats, nT);
        E     = reshape(ve(1,:,:), num_sats, nT);
        N     = reshape(ve(2,:,:), num_sats, nT);
        U     = reshape(ve(3,:,:), num_sats, nT);
        r_m   = sqrt(E.^2 + N.^2 + U.^2);
        el_m  = asind(U ./ r_m);
        r_m(el_m < min_elev) = Inf;
        [best_r, bi] = min(r_m, [], 1);
        has_srv   = ~isinf(best_r);
        best_el_i  = NaN(1, nT);
        best_rng_i = NaN(1, nT);
        best_rng_i(has_srv) = best_r(has_srv);
        vc = find(has_srv);
        best_el_i(has_srv) = el_m(bi(has_srv) + (vc - 1) * num_sats);
        best_el(i, :)  = best_el_i;
        best_rng(i, :) = best_rng_i;
    end
end

function run_aer(ue_gs, sats, min_elev)
% Native toolbox aer(), per-UE loop.
    num_ues  = numel(ue_gs);
    [~, el0, ~] = aer(ue_gs(1), sats);
    nT = size(el0, 2);
    nS = size(el0, 1);
    best_el  = NaN(num_ues, nT);
    best_rng = NaN(num_ues, nT);
    for i = 1:num_ues
        [~, el_mat, r_mat] = aer(ue_gs(i), sats);
        r_mat(el_mat < min_elev) = Inf;
        [best_r, bi] = min(r_mat, [], 1);
        has_srv = ~isinf(best_r);
        best_rng(i, has_srv) = best_r(has_srv);
        vc = find(has_srv);
        best_el(i, has_srv) = el_mat(bi(has_srv) + (vc - 1) * nS);
    end
end

function run_vec(ue_lats, ue_lons, sat_ecef, min_elev)
% Vectorised ECEF geometry, sequential per-UE loop.
    ue_xyz   = lla2ecef([ue_lats(:), ue_lons(:), zeros(numel(ue_lats), 1)]);
    num_ues  = numel(ue_lats);
    num_sats = size(sat_ecef, 2);
    nT       = size(sat_ecef, 3);
    best_el  = NaN(num_ues, nT);
    best_rng = NaN(num_ues, nT);
    for i = 1:num_ues
        dx    = sat_ecef - ue_xyz(i, :)';
        slat  = sind(ue_lats(i));  clat = cosd(ue_lats(i));
        slon  = sind(ue_lons(i));  clon = cosd(ue_lons(i));
        R_enu = [-slon,       clon,       0; ...
                 -slat*clon, -slat*slon,  clat; ...
                  clat*clon,  clat*slon,  slat];
        ve    = reshape(R_enu * reshape(dx, 3, []), 3, num_sats, nT);
        E     = reshape(ve(1,:,:), num_sats, nT);
        N     = reshape(ve(2,:,:), num_sats, nT);
        U     = reshape(ve(3,:,:), num_sats, nT);
        r_m   = sqrt(E.^2 + N.^2 + U.^2);
        el_m  = asind(U ./ r_m);
        r_m(el_m < min_elev) = Inf;
        [best_r, bi] = min(r_m, [], 1);
        has_srv = ~isinf(best_r);
        best_rng(i, has_srv) = best_r(has_srv);
        vc = find(has_srv);
        best_el(i, has_srv) = el_m(bi(has_srv) + (vc - 1) * num_sats);
    end
end


