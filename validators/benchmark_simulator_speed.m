% benchmark_simulator_speed.m
% Combined benchmark: AER geometry methods, parallelism case studies,
% and Constellation_simulator worker scaling.
%
% Section A: AER method comparison — geometry cost vs UE count
%     1. aer()          — per-UE toolbox call (baseline)
%     2. vec-toolbox    — states() ECEF prefetch + vectorised ENU loop
%     3. vec-fastmath   — fast_walker_ecef() + same ENU loop
%
% Section B: Parallelism case studies
%     4. par-gridsearch — many independent small sims via parfor
%     5. par-toolbox    — parfor over UEs within a single large sim
%
% Section C: Constellation_simulator worker scaling
%     Toolbox propagator (Processes pool, matches Stage-3 gridsearch).
%     Sweeps parfeval workers 1→32. Shows throughput (sims/s) vs simultaneous workers.
%     Parameters match gridsearch Stage-3 defaults.
%
% Figures saved to: plotting_scripts/figures/optimisations/

clear; close all; clc;
% Add the project to the MATLAB path (robust to the script's folder depth).
repo_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
addpath(fullfile(repo_root, 'functions'));
path_setup();

out_dir = fullfile(repo_root, 'plotting_scripts', 'figures', 'optimisations');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% =========================================================================
% SHARED CONSTELLATION PARAMETERS
% =========================================================================
height_km = 1000;
height_m  = height_km * 1e3;
r_earth_m = 6378.137e3;
Min_elev  = 20;
Lat_range = [54 + 35/60, 83 + 40/60];   % Denmark–Greenland

opt     = load(fullfile(fileparts(mfilename('fullpath')), 'optimal_constellations.mat'));
[~, ci] = min(abs(opt.heights_km - height_km));
con     = table2struct(opt.best_delta_sats(ci, :));

Total_sats     = con.Sats_per_plane * con.Num_planes;
Num_planes     = con.Num_planes;
Sats_per_plane = con.Sats_per_plane;
Phasing        = con.Phasing;
Inclination    = con.Inclination;

start_time = datetime(2025, 6, 1, 12, 0, 0, 'TimeZone', 'UTC');

fprintf('=== Constellation Simulator Speed Benchmark ===\n');
fprintf('Constellation: %d sats | %d planes | %.0f deg inc | %d km\n\n', ...
    Total_sats, Num_planes, Inclination, height_km);

% =========================================================================
% SECTION A: AER method comparison — geometry cost vs UE count
%   Fixed 2500 timesteps at 60 s.  Sweeps three geometry methods over
%   increasing UE counts.
% =========================================================================
fprintf('--- Section A: AER method comparison ---\n');

n_steps_aer = 2500;
sample_aer  = 60;   % seconds
stop_aer    = start_time + seconds((n_steps_aer - 1) * sample_aer);
t_steps_aer = 0 : sample_aer : (n_steps_aer - 1) * sample_aer;

ue_counts = [50, 150, 500, 1000, 2500];
n_tests   = numel(ue_counts);

rng(42);
n_side = ceil(sqrt(max(ue_counts) * 2));
[LAT, LON] = meshgrid(linspace(Lat_range(1), Lat_range(2), n_side), ...
                       linspace(-60, 30, n_side));
pool_lats = LAT(:);  pool_lons = LON(:);
perm_idx  = randperm(numel(pool_lats));

fprintf('  Pre-computing sat positions (fast_walker_ecef)... ');
tw = tic;
sat_ecef_fast = fast_walker_ecef(height_m, Inclination, Num_planes, ...
    Sats_per_plane, Phasing, t_steps_aer, start_time);
t_fast_setup = toc(tw);
fprintf('%.2f s  (%d sats x %d steps)\n\n', ...
    t_fast_setup, size(sat_ecef_fast,2), size(sat_ecef_fast,3));

t_setup = zeros(3, n_tests);
t_geom  = zeros(3, n_tests);
t_total = zeros(3, n_tests);

for ti = 1:n_tests
    n_ues   = ue_counts(ti);
    ue_lats = pool_lats(perm_idx(1:n_ues));
    ue_lons = pool_lons(perm_idx(1:n_ues));

    fprintf('  [%d/%d]  %4d UEs x %d steps\n', ti, n_tests, n_ues, n_steps_aer);

    % Method 1: native aer()
    tw = tic;  ts = tic;
        sc1 = satelliteScenario;
        sc1.StartTime = start_time;  sc1.StopTime = stop_aer;  sc1.SampleTime = sample_aer;
        sats1 = walkerDelta(sc1, height_m + r_earth_m, Inclination, ...
                            Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
        ue_gs = groundStation(sc1, ue_lats, ue_lons);
    t_setup(1,ti) = toc(ts);  ts = tic;
        run_aer(ue_gs, sats1, Min_elev);
    t_geom(1,ti)  = toc(ts);
    t_total(1,ti) = toc(tw);
    delete(sc1);

    % Method 2: vec-toolbox
    tw = tic;  ts = tic;
        sc2 = satelliteScenario;
        sc2.StartTime = start_time;  sc2.StopTime = stop_aer;  sc2.SampleTime = sample_aer;
        sats2 = walkerDelta(sc2, height_m + r_earth_m, Inclination, ...
                            Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
        [raw, ~, ~] = states(sats2, "CoordinateFrame", "ECEF");
        sat_ecef_tb = permute(raw, [1, 3, 2]);
    t_setup(2,ti) = toc(ts);  ts = tic;
        run_vec(ue_lats, ue_lons, sat_ecef_tb, Min_elev);
    t_geom(2,ti)  = toc(ts);
    t_total(2,ti) = toc(tw);
    delete(sc2);

    % Method 3: vec-fastmath (sat positions pre-computed above, reused)
    t_setup(3,ti) = t_fast_setup;
    ts = tic;
        run_vec(ue_lats, ue_lons, sat_ecef_fast, Min_elev);
    t_geom(3,ti)  = toc(ts);
    t_total(3,ti) = t_fast_setup + t_geom(3,ti);

    fprintf('    aer():       %6.2f s  (setup %5.2f + geom %5.2f)\n', t_total(1,ti), t_setup(1,ti), t_geom(1,ti));
    fprintf('    vec-toolbox: %6.2f s  (setup %5.2f + geom %5.2f)\n', t_total(2,ti), t_setup(2,ti), t_geom(2,ti));
    fprintf('    vec-fastmath:%6.2f s  (setup %5.2f + geom %5.2f)\n\n', t_total(3,ti), t_fast_setup, t_geom(3,ti));
end

% =========================================================================
% SECTION B: Parallelism case studies
%   4. par-gridsearch: n_sims independent small sims via parfor
%      (models real gridsearch: many configs × few UEs each)
%   5. par-toolbox:    parfor over UEs within a single large sim
%      (models Constellation_simulator with use_parallel=true)
% =========================================================================
fprintf('--- Section B: Parallelism case studies ---\n');

% Start pool early (needed for both B and C)
pool = gcp('nocreate');
if isempty(pool), pool = parpool('Threads'); end

% Method 4: par-gridsearch
n_ues_cs     = 50;
n_per_worker = 10;
n_sims       = n_per_worker * 32;   % = 320

ue_lats_cs = cell(n_sims, 1);
ue_lons_cs = cell(n_sims, 1);
for k = 1:n_sims
    i0 = mod((k-1) * n_ues_cs, numel(pool_lats) - n_ues_cs) + 1;
    ue_lats_cs{k} = pool_lats(perm_idx(i0 : i0 + n_ues_cs - 1));
    ue_lons_cs{k} = pool_lons(perm_idx(i0 : i0 + n_ues_cs - 1));
end

fprintf('  Method 4 — sequential baseline: %d x %d-UE fastmath sims...\n', n_sims, n_ues_cs);
ts = tic;
for k = 1:n_sims
    run_vec(ue_lats_cs{k}, ue_lons_cs{k}, sat_ecef_fast, Min_elev);
end
t_seq_total = toc(ts);
fprintf('    sequential: %.2f s  (%.3f s/sim)\n\n', t_seq_total, t_seq_total/n_sims);

% Warm-up parfor
parfor k = 1:2
    run_vec(ue_lats_cs{k}, ue_lons_cs{k}, sat_ecef_fast, Min_elev);
end

fprintf('  Method 4 — parallel: %d x %d-UE fastmath sims (parfor, ~%d/worker)...\n', ...
    n_sims, n_ues_cs, n_per_worker);
ts = tic;
parfor k = 1:n_sims
    run_vec(ue_lats_cs{k}, ue_lons_cs{k}, sat_ecef_fast, Min_elev);
end
t_par_gridsearch = toc(ts);
par_speedup = t_seq_total / t_par_gridsearch;
fprintf('    parallel: %.2f s  (%.3f s/sim)  ->  speedup: %.1fx\n\n', ...
    t_par_gridsearch, t_par_gridsearch/n_sims, par_speedup);

% Method 5: par-toolbox (parfor over UEs, single large sim)
n_ues_par5   = 2500;
ue_lats_par5 = pool_lats(perm_idx(1:n_ues_par5));
ue_lons_par5 = pool_lons(perm_idx(1:n_ues_par5));

fprintf('  Method 5 — par-toolbox: %d UEs x %d steps...\n', n_ues_par5, n_steps_aer);
tw = tic;  ts = tic;
    sc5 = satelliteScenario;
    sc5.StartTime = start_time;  sc5.StopTime = stop_aer;  sc5.SampleTime = sample_aer;
    sats5 = walkerDelta(sc5, height_m + r_earth_m, Inclination, ...
                        Total_sats, Num_planes, Phasing, OrbitPropagator="sgp4");
    [raw5, ~, ~]    = states(sats5, "CoordinateFrame", "ECEF");
    sat_ecef_par_tb = permute(raw5, [1, 3, 2]);
t_par_tb_setup = toc(ts);  ts = tic;
    run_vec_par_ues(ue_lats_par5, ue_lons_par5, sat_ecef_par_tb, Min_elev);
t_par_tb_geom  = toc(ts);
t_par_tb_total = toc(tw);
delete(sc5);

tb_seq_total = t_setup(2,end) + t_geom(2,end);
tb_speedup   = tb_seq_total / t_par_tb_total;
fprintf('    par-toolbox: %.2f s  (setup %.2f + geom %.2f)  ->  speedup vs vec-toolbox: %.1fx\n\n', ...
    t_par_tb_total, t_par_tb_setup, t_par_tb_geom, tb_speedup);

% =========================================================================
% SECTION C: Constellation_simulator worker scaling
%   Vary simultaneous parfeval calls 1→32.  Toolbox propagator (use_toolbox=true),
%   matching Stage-3 gridsearch.  Requires a Processes pool — toolbox objects
%   (satelliteScenario) cannot be constructed on Threads pool workers.
%     SampleTime=660 s,  nUEs=nSteps≈2146  (certainty=99%, fa=ft=0.1%)
% =========================================================================
fprintf('--- Section C: Constellation_simulator worker scaling ---\n');

certainty    = 0.99;
fa = 0.001;  ft = 0.001;
n_ues_stage3 = ceil(sqrt(log(1 - certainty) / log(1 - fa * ft)));   % ≈ 2146
sample_time  = 660;   % seconds
n_steps      = n_ues_stage3;
stop_stage3  = start_time + seconds((n_steps - 1) * sample_time);

fprintf('  Stage-3: %d UEs x %d steps x %d sats\n\n', n_ues_stage3, n_steps, Total_sats);

[lats_c, lons_c] = generate_equal_area_ues(Lat_range, [-180, 180], n_ues_stage3);
Cfg.StartTime          = start_time;
Cfg.StopTime           = stop_stage3;
Cfg.SampleTime         = sample_time;
Cfg.Orbit_height       = height_m;
Cfg.Inclination        = Inclination;
Cfg.Num_planes         = Num_planes;
Cfg.Sats_per_plane     = Sats_per_plane;
Cfg.Total_sats         = Total_sats;
Cfg.Phasing            = Phasing;
Cfg.Min_elevation_UE   = Min_elev;
Cfg.Lat_range_deg      = Lat_range;
Cfg.WalkerStar         = false;
Cfg.Flat_UE_array.Lats = lats_c;
Cfg.Flat_UE_array.Lons = lons_c;

% Section C requires a Processes pool so each worker can instantiate
% satelliteScenario objects.  Switch away from the Threads pool used in B.
if ~isempty(pool) && isa(pool, 'parallel.ThreadPool')
    fprintf('  Switching from Threads to Processes pool for Section C...\n');
    delete(pool);
    pool = parpool();
end

N_WORKERS = 32;
nw        = pool.NumWorkers;

% Warm up pool
fprintf('  Warming up pool (%d workers)...\n', nw);
wf = cell(nw, 1);
for w = 1:nw
    wf{w} = parfeval(pool, @(C) Constellation_simulator(C, false, false, true), 1, Cfg);
end
for w = 1:nw
    wait(wf{w});
    if ~isempty(wf{w}.Error)
        warning('Warm-up worker %d failed: %s', w, wf{w}.Error.message);
    end
end
fprintf('  Warm-up done.\n\n');

worker_counts = unique([1, 2, 4, 8, 16, N_WORKERS]);
worker_counts = worker_counts(worker_counts <= nw);
n_wc = numel(worker_counts);

t_toolbox_c = zeros(1, n_wc);

for ci = 1:n_wc
    npar = worker_counts(ci);

    futs = cell(npar, 1);  tw = tic;
    for w = 1:npar
        futs{w} = parfeval(pool, @(C) Constellation_simulator(C, false, false, true), 1, Cfg);
    end
    for w = 1:npar
        wait(futs{w});
        if ~isempty(futs{w}.Error)
            warning('Worker %d failed: %s', w, futs{w}.Error.message);
        end
    end
    t_toolbox_c(ci) = toc(tw);

    fprintf('  %2d workers:  %.2f s  (%.3f sims/s)\n', npar, t_toolbox_c(ci), npar/t_toolbox_c(ci));
end

thpt_toolbox   = worker_counts ./ t_toolbox_c;
[pk_tb, pi_tb] = max(thpt_toolbox);
fprintf('\n  Peak throughput: %.2f sims/s at %d workers\n\n', pk_tb, worker_counts(pi_tb));

% =========================================================================
% SUMMARY TABLES
% =========================================================================
speedup_tb   = t_total(1,:) ./ t_total(2,:);
speedup_fast = t_total(1,:) ./ t_total(3,:);

W = 78;
fprintf('%s\n', repmat('=', 1, W));
fprintf('  %-14s  %6s  %10s  %10s  %10s  %9s\n', 'Method', 'UEs', 'Setup (s)', 'Geom (s)', 'Total (s)', 'Speedup');
fprintf('%s\n', repmat('-', 1, W));
names = {'aer()', 'vec-toolbox', 'vec-fastmath'};
for ti = 1:n_tests
    for m = 1:3
        fprintf('  %-14s  %6d  %10.2f  %10.2f  %10.2f  %8.1fx\n', ...
            names{m}, ue_counts(ti), t_setup(m,ti), t_geom(m,ti), t_total(m,ti), t_total(1,ti)/t_total(m,ti));
    end
    if ti < n_tests, fprintf('%s\n', repmat('-', 1, W)); end
end
fprintf('%s\n', repmat('=', 1, W));
fprintf('\n  par-toolbox case study: %d UEs x %d steps\n', n_ues_par5, n_steps_aer);
fprintf('  vec-toolbox (seq): %.2f s  |  par-toolbox: %.2f s  |  Speedup: %.1fx\n', ...
    tb_seq_total, t_par_tb_total, tb_speedup);
fprintf('\n  Gridsearch case study: %d independent %d-UE fastmath sims\n', n_sims, n_ues_cs);
fprintf('  Sequential: %.2f s  |  Parallel: %.2f s  |  Speedup: %.1fx\n', ...
    t_seq_total, t_par_gridsearch, par_speedup);
fprintf('%s\n\n', repmat('=', 1, W));

fprintf('%s\n', repmat('=', 1, 52));
fprintf('  %8s  %8s  %10s\n', 'Workers', 'Time (s)', 'Thpt (sims/s)');
fprintf('%s\n', repmat('-', 1, 52));
for ci = 1:n_wc
    fprintf('  %8d  %8.2f  %10.3f\n', worker_counts(ci), t_toolbox_c(ci), thpt_toolbox(ci));
end
fprintf('%s\n\n', repmat('=', 1, 52));

% =========================================================================
% FIGURES
% =========================================================================
c1 = [0.850, 0.325, 0.098];
c2 = [0,     0.447, 0.741];
c3 = [0.466, 0.674, 0.188];

% --- Figure 1: Computation time vs UE count (log-log) --------------------
f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 620, 420]);
ax1 = axes(f1, 'Color', 'w');  hold(ax1, 'on');

semilogy(ax1, ue_counts, t_total(1,:), '-o',  'LineWidth', 2.5, 'Color', c1,     'MarkerFaceColor', c1,     'MarkerSize', 7, 'DisplayName', 'aer()  — total');
semilogy(ax1, ue_counts, t_geom(1,:),  '--o', 'LineWidth', 1.5, 'Color', c1*0.6, 'MarkerSize', 5,           'DisplayName', 'aer()  — geom only');
semilogy(ax1, ue_counts, t_total(2,:), '-s',  'LineWidth', 2.5, 'Color', c2,     'MarkerFaceColor', c2,     'MarkerSize', 7, 'DisplayName', 'vec-toolbox  — total');
semilogy(ax1, ue_counts, t_geom(2,:),  '--s', 'LineWidth', 1.5, 'Color', c2*0.6, 'MarkerSize', 5,           'DisplayName', 'vec-toolbox  — geom only');
semilogy(ax1, ue_counts, t_total(3,:), '-^',  'LineWidth', 2.5, 'Color', c3,     'MarkerFaceColor', c3,     'MarkerSize', 7, 'DisplayName', 'vec-fastmath  — total');
semilogy(ax1, ue_counts, t_geom(3,:),  '--^', 'LineWidth', 1.5, 'Color', c3*0.6, 'MarkerSize', 5,           'DisplayName', 'vec-fastmath  — geom only');

ax1.XScale = 'log';
xlabel(ax1, 'Number of UEs',        'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax1, 'Computation time (s)', 'FontWeight', 'bold', 'FontSize', 13);
legend(ax1, 'Location', 'northwest', 'FontSize', 9);
grid(ax1, 'on');  ax1.Box = 'off';
title(ax1, sprintf('%d sats @ %d km  |  %d timesteps  (\\DeltaT = %d s)', ...
    Total_sats, height_km, n_steps_aer, sample_aer), 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f1, fullfile(out_dir, 'aer_speedup_timing.png'), 'Resolution', 300);
close(f1);  fprintf('Saved: aer_speedup_timing.png\n');

% --- Figure 2: Speedup factor vs UE count --------------------------------
f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 500, 360]);
ax2 = axes(f2, 'Color', 'w');  hold(ax2, 'on');

plot(ax2, ue_counts, speedup_tb,   '-s', 'LineWidth', 2.5, 'Color', c2, 'MarkerFaceColor', c2, 'MarkerSize', 8, 'DisplayName', 'vec-toolbox');
plot(ax2, ue_counts, speedup_fast, '-^', 'LineWidth', 2.5, 'Color', c3, 'MarkerFaceColor', c3, 'MarkerSize', 8, 'DisplayName', 'vec-fastmath');
yline(ax2, 1, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

ax2.XScale = 'log';
xlabel(ax2, 'Number of UEs',                   'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax2, 'Speedup  (t_{aer} / t_{method})', 'FontWeight', 'bold', 'FontSize', 13);
legend(ax2, 'Location', 'northwest', 'FontSize', 11);
grid(ax2, 'on');  ax2.Box = 'off';
title(ax2, sprintf('Speedup relative to aer()  |  %d sats, %d timesteps', Total_sats, n_steps_aer), ...
    'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f2, fullfile(out_dir, 'aer_speedup_factor.png'), 'Resolution', 300);
close(f2);  fprintf('Saved: aer_speedup_factor.png\n');

% --- Figure 3: Parallelism case study (stacked bar) ----------------------
bar_setup4 = [t_setup(2,end), t_par_tb_setup,  t_fast_setup,      t_fast_setup];
bar_geom4  = [t_geom(2,end),  t_par_tb_geom,   t_seq_total,       t_par_gridsearch];
bar_labels = {'vec-toolbox', 'par-toolbox', ...
              sprintf('seq (%d\\times%d UE)', n_sims, n_ues_cs), ...
              sprintf('par (%d\\times%d UE)', n_sims, n_ues_cs)};

f3 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 600, 390]);
ax3 = axes(f3, 'Color', 'w');

bh = bar(ax3, [bar_setup4; bar_geom4]', 'stacked');
bh(1).FaceColor = [0.72, 0.72, 0.72];  bh(1).DisplayName = 'Setup / orbit generation';
bh(2).FaceColor = [0.18, 0.38, 0.60];  bh(2).DisplayName = 'Geometry (total)';

totals4  = bar_setup4 + bar_geom4;
max_bar4 = max(totals4);
text(ax3, 1, totals4(1) + max_bar4*0.03, '1\times',                           'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
text(ax3, 2, totals4(2) + max_bar4*0.03, sprintf('%.0f\\times', tb_speedup),  'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', c2);
text(ax3, 3, totals4(3) + max_bar4*0.03, '1\times',                           'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
text(ax3, 4, totals4(4) + max_bar4*0.03, sprintf('%.0f\\times', par_speedup), 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', c3);

xline(ax3, 2.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
ax3.XTickLabel = bar_labels;  ax3.XTickLabelRotation = 8;
ylabel(ax3, 'Total wall-clock time (s)', 'FontWeight', 'bold', 'FontSize', 13);
legend(ax3, 'Location', 'northeast', 'FontSize', 11);
grid(ax3, 'on');  ax3.Box = 'off';  ax3.XGrid = 'off';
title(ax3, sprintf('%d sats @ %d km  |  left: 2500 UEs\\times%d steps  |  right: %d\\times%d-UE sims', ...
    Total_sats, height_km, n_steps_aer, n_sims, n_ues_cs), 'FontSize', 10, 'FontWeight', 'bold');

exportgraphics(f3, fullfile(out_dir, 'aer_speedup_case_study.png'), 'Resolution', 300);
close(f3);  fprintf('Saved: aer_speedup_case_study.png\n');

% --- Figure 4: Wall-clock time vs simultaneous workers -------------------
f4 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 680, 440]);
ax4 = axes(f4, 'Color', 'w');  hold(ax4, 'on');

plot(ax4, worker_counts, t_toolbox_c, '-s', 'LineWidth', 2.5, 'Color', c2, 'MarkerFaceColor', c2, 'MarkerSize', 8, 'DisplayName', 'toolbox propagator');

x_ref = linspace(worker_counts(1), worker_counts(end), 200);
plot(ax4, x_ref, t_toolbox_c(1) .* x_ref, ':', 'LineWidth', 1.2, 'Color', c2*0.65, 'DisplayName', 'linear (no speedup)');

xlabel(ax4, 'Simultaneous workers', 'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax4, 'Wall-clock time (s)',  'FontWeight', 'bold', 'FontSize', 13);
legend(ax4, 'Location', 'northwest', 'FontSize', 11);
grid(ax4, 'on');  ax4.Box = 'off';
title(ax4, sprintf('Stage-3 toolbox: %d UEs\\times%d steps\\times%d sats', ...
    n_ues_stage3, n_steps, Total_sats), 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f4, fullfile(out_dir, 'simulator_wall_time_vs_workers.png'), 'Resolution', 300);
close(f4);  fprintf('Saved: simulator_wall_time_vs_workers.png\n');

% --- Figure 5: Throughput (sims/s) vs simultaneous workers ---------------
f5 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 680, 440]);
ax5 = axes(f5, 'Color', 'w');  hold(ax5, 'on');

plot(ax5, worker_counts, thpt_toolbox, '-s', 'LineWidth', 2.5, 'Color', c2, 'MarkerFaceColor', c2, 'MarkerSize', 8, 'DisplayName', 'toolbox propagator');
plot(ax5, worker_counts, worker_counts / t_toolbox_c(1), '--', 'LineWidth', 1.2, 'Color', c2*0.65, 'DisplayName', 'linear ideal');

text(ax5, worker_counts(pi_tb), pk_tb * 1.06, sprintf('%.2f sims/s', pk_tb), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', c2);

xlabel(ax5, 'Simultaneous workers',       'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax5, 'Throughput (sims / second)', 'FontWeight', 'bold', 'FontSize', 13);
legend(ax5, 'Location', 'northwest',      'FontSize', 11);
grid(ax5, 'on');  ax5.Box = 'off';
title(ax5, sprintf('Gridsearch throughput  |  %d UEs\\times%d steps  |  %d sats', ...
    n_ues_stage3, n_steps, Total_sats), 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f5, fullfile(out_dir, 'simulator_throughput_vs_workers.png'), 'Resolution', 300);
close(f5);  fprintf('Saved: simulator_throughput_vs_workers.png\n');

fprintf('\nAll figures saved to: %s\n', out_dir);

% =========================================================================
% LATEX TABLE: Worker scaling
% =========================================================================
tex_workers = fullfile(out_dir, 'worker_scaling_table.tex');
fid = fopen(tex_workers, 'w');
fprintf(fid, '\\begin{table}[H]\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{\\texttt{constellation\\_simulator} wall-clock time and throughput as a function\n');
fprintf(fid, '             of simultaneously dispatched workers (\\texttt{parfeval}, Processes pool).\n');
fprintf(fid, '             Stage-3 parameters: %d \\acp{UE} $\\times$ %d timesteps $\\times$ %d satellites.\n', ...
    n_ues_stage3, n_steps, Total_sats);
fprintf(fid, '             Toolbox propagator (two-body-Keplerian).}\n');
fprintf(fid, '    \\label{tab:worker_scaling}\n');
fprintf(fid, '    \\begin{tabular}{r r r r}\n');
fprintf(fid, '        \\toprule\n');
fprintf(fid, '        \\textbf{Workers} & \\textbf{Wall-clock (s)} & \\textbf{Throughput (sims/s)} & \\textbf{Efficiency} \\\\\n');
fprintf(fid, '        \\midrule\n');
for ci = 1:n_wc
    efficiency = (thpt_toolbox(ci) / (worker_counts(ci) / t_toolbox_c(1))) * 100;
    fprintf(fid, '        %2d & %.2f & %.3f & %.0f\\%% \\\\\n', ...
        worker_counts(ci), t_toolbox_c(ci), thpt_toolbox(ci), efficiency);
end
fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('Saved: worker_scaling_table.tex\n');

% =========================================================================
% LATEX TABLE: AER method speedup
% =========================================================================
latex_path = fullfile(out_dir, 'aer_speedup_table.tex');
fid = fopen(latex_path, 'w');
fprintf(fid, '\\begin{table}[H]\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Azimuth Elevation and Range, AER, geometry computation time as a function of \\ac{UE} count.\n');
fprintf(fid, '             Constellation: %d satellites at %d\\,km, %d timesteps.\n', ...
    Total_sats, height_km, n_steps_aer);
fprintf(fid, '             Speedup is relative to the MATLAB \\texttt{aer()} baseline.}\n');
fprintf(fid, '    \\label{tab:aer_speedup}\n');
fprintf(fid, '    \\begin{tabular}{r r r r r r}\n');
fprintf(fid, '        \\toprule\n');
fprintf(fid, '        & \\multicolumn{1}{c}{\\textbf{aer()}}\n');
fprintf(fid, '        & \\multicolumn{2}{c}{\\textbf{Custom AER}}\n');
fprintf(fid, '        & \\multicolumn{2}{c}{\\textbf{+ custom propagation}} \\\\\n');
fprintf(fid, '        \\cmidrule(lr){2-2}\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n');
fprintf(fid, '        \\textbf{UEs} & \\textbf{Total (s)} & \\textbf{Total (s)} & \\textbf{Speedup}\n');
fprintf(fid, '                     & \\textbf{Total (s)} & \\textbf{Speedup} \\\\\n');
fprintf(fid, '        \\midrule\n');
for ti = 1:n_tests
    su_tb   = round(t_total(1,ti) / t_total(2,ti));
    su_fast = round(t_total(1,ti) / t_total(3,ti));
    fprintf(fid, '        %5d & %6.2f & %5.2f & $%d\\times$ & %5.2f & $%d\\times$ \\\\\n', ...
        ue_counts(ti), t_total(1,ti), t_total(2,ti), su_tb, t_total(3,ti), su_fast);
end
fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('Saved: aer_speedup_table.tex\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function run_aer(ue_gs, sats, min_elev)
% Native toolbox aer(), per-UE loop.
    num_ues = numel(ue_gs);
    [~, el0, ~] = aer(ue_gs(1), sats);
    nT = size(el0, 2);  nS = size(el0, 1);
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

function run_vec_par_ues(ue_lats, ue_lons, sat_ecef, min_elev)
% Vectorised ECEF geometry, parfor over UEs.
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
        has_srv = ~isinf(best_r);
        best_el_i  = NaN(1, nT);
        best_rng_i = NaN(1, nT);
        best_rng_i(has_srv) = best_r(has_srv);
        vc = find(has_srv);
        best_el_i(has_srv) = el_m(bi(has_srv) + (vc - 1) * num_sats);
        best_el(i, :)  = best_el_i;
        best_rng(i, :) = best_rng_i;
    end
end
