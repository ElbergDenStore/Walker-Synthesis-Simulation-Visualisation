% compare_simulator_speed.m
% Measures and explains the performance gap between:
%
%   1. coverage_simulator_function          (original)  — struct-per-UE writes, double
%   2. fast_coverage_simulator_function_v2  (v2)        — plain-matrix writes, double
%   3. fast_coverage_simulator_function_v3  (v3)        — plain-matrix writes, single
%                                             sat_pos_ecef cast to single so 1.44 MB
%                                             fits in P-core L2 (2 MB) instead of
%                                             spilling to L3 every UE iteration.
%
% Three scenarios that match real gridsearch conditions:
%
%   A. Single worker       — one call, isolates struct overhead from cache effects
%   B. 32 parallel workers — parfeval × N_WORKERS simultaneous calls (gridsearch)
%   C. Worker sweep 1→32   — shows L3 cache thrash knee per function

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions', 'data'));

out_dir = fullfile(fileparts(mfilename('fullpath')), 'plotting_scripts/figures', 'optimisations');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% =========================================================================
% SIMULATION PARAMETERS  — mirror gridsearch Stage-3 defaults
%   (run_single_gridsearch: SampleTime=660, certainty=99%, fa=0.1%, ft=0.1%)
% =========================================================================
N_WORKERS = 32;                          % workers to launch for parallel scenario

certainty = 0.99;
fa = 0.001;  ft = 0.001;
req_samples = log(1-certainty) / log(1-fa*ft);
n_ues_stage3 = ceil(sqrt(req_samples));                 % = 2146
sample_time  = 660;                                     % seconds
n_steps      = n_ues_stage3;                            % square: nUEs = nSteps = 2146

height_m     = 1000e3;
Inclination  = 77;
Num_planes   = 4;
Sats_per_plane = 14;
Phasing      = 2;
Lat_range    = [54+35/60, 83+40/60];   % Denmark–Greenland
Min_elev     = 20;

fprintf('=== Simulator Speed Comparison ===\n');
fprintf('Stage-3 equivalent: %d UEs x %d steps x %d sats\n\n', ...
    n_ues_stage3, n_steps, Num_planes*Sats_per_plane);

% =========================================================================
% BUILD SHARED Cfg (used by both functions)
% =========================================================================
start_time = datetime('2025-06-01 12:00:00', 'TimeZone', 'UTC');
stop_time  = start_time + seconds((n_steps-1) * sample_time);

[lats, lons] = generate_equal_ish_area_UEs(Lat_range, [-180, 180], n_ues_stage3);

Cfg.StartTime        = start_time;
Cfg.StopTime         = stop_time;
Cfg.SampleTime       = sample_time;
Cfg.Orbit_height     = height_m;
Cfg.Inclination      = Inclination;
Cfg.Num_planes       = Num_planes;
Cfg.Sats_per_plane   = Sats_per_plane;
Cfg.Total_sats       = Num_planes * Sats_per_plane;
Cfg.Phasing          = Phasing;
Cfg.Min_elevation_UE = Min_elev;
Cfg.Lat_range_deg    = Lat_range;
Cfg.WalkerStar       = false;
Cfg.Flat_UE_array.Lats = lats;
Cfg.Flat_UE_array.Lons = lons;

% =========================================================================
% SCENARIO A: Single worker — one call at a time
%   Measures raw single-call latency without any cache or lock contention.
%   N_REP repetitions to get a stable average.
% =========================================================================
N_REP = 3;
fprintf('--- Scenario A: single worker, %d repetitions ---\n', N_REP);

t_old_single = zeros(N_REP, 1);
t_new_single = zeros(N_REP, 1);
t_v3_single  = zeros(N_REP, 1);

t_v3_single = zeros(N_REP, 1);

for r = 1:N_REP
    t0 = tic;
    coverage_simulator_function(Cfg, false, false, []);
    t_old_single(r) = toc(t0);
    fprintf('  old v1 rep %d: %.2f s\n', r, t_old_single(r));
end
for r = 1:N_REP
    t0 = tic;
    fast_coverage_simulator_function_v2(Cfg, false, false, false, []);
    t_new_single(r) = toc(t0);
    fprintf('  new v2 rep %d: %.2f s\n', r, t_new_single(r));
end
for r = 1:N_REP
    t0 = tic;
    fast_coverage_simulator_function_v3(Cfg, false, false, false, []);
    t_v3_single(r) = toc(t0);
    fprintf('  new v3 rep %d: %.2f s\n', r, t_v3_single(r));
end

t_old_A = median(t_old_single);
t_new_A = median(t_new_single);
t_v3_A  = median(t_v3_single);
speedup_A_v2 = t_old_A / t_new_A;
speedup_A_v3 = t_old_A / t_v3_A;
fprintf('  >> v2 speedup: %.1fx  (old %.2fs  vs  v2 %.2fs)\n',   speedup_A_v2, t_old_A, t_new_A);
fprintf('  >> v3 speedup: %.1fx  (old %.2fs  vs  v3 %.2fs)\n\n', speedup_A_v3, t_old_A, t_v3_A);

% =========================================================================
% SCENARIO B: 32 parallel workers — parfeval × N_WORKERS simultaneous calls
%   This is exactly what gridsearch does for Stage 3.
%   Each worker is a separate MATLAB process; they share the L3 cache and
%   memory bus but have independent heaps.
%
%   Timing: wall clock from first parfeval dispatch to last fetchOutputs.
% =========================================================================
fprintf('--- Scenario B: %d simultaneous workers (gridsearch-equivalent) ---\n', N_WORKERS);

pool = gcp('nocreate');
if isempty(pool), pool = parpool(); end
nw = pool.NumWorkers;
fprintf('  Pool has %d workers (requested %d).\n', nw, N_WORKERS);

% Warm up workers (one dummy call each) to ensure any one-time startup cost
% is excluded from the benchmark timing.
fprintf('  Warming up workers...\n');
wf = cell(nw, 1);
for w = 1:nw
    wf{w} = parfeval(pool, @(C) coverage_simulator_function(C, false, false, []), 0, Cfg);
end
for w = 1:nw, wait(wf{w}); end
fprintf('  Warm-up complete.\n\n');

% -- Old v1, N_WORKERS simultaneous --
fprintf('  Dispatching %d × old v1...\n', N_WORKERS);
t_wall_old_B = tic;
futs_old = cell(N_WORKERS, 1);
for w = 1:N_WORKERS
    futs_old{w} = parfeval(pool, @(C) coverage_simulator_function(C, false, false, []), 0, Cfg);
end
for w = 1:N_WORKERS, wait(futs_old{w}); end
t_old_B = toc(t_wall_old_B);
fprintf('  old v1: %.2f s wall  (%.2f s per worker)\n', t_old_B, t_old_B);

% -- New v2, N_WORKERS simultaneous --
fprintf('  Dispatching %d × new v2...\n', N_WORKERS);
t_wall_new_B = tic;
futs_new = cell(N_WORKERS, 1);
for w = 1:N_WORKERS
    futs_new{w} = parfeval(pool, @(C) fast_coverage_simulator_function_v2(C, false, false, false, []), 0, Cfg);
end
for w = 1:N_WORKERS, wait(futs_new{w}); end
t_new_B = toc(t_wall_new_B);
fprintf('  new v2: %.2f s wall  (%.2f s per worker)\n', t_new_B, t_new_B);

% -- New v3 (single-precision), N_WORKERS simultaneous --
fprintf('  Dispatching %d × new v3 (single)...\n', N_WORKERS);
t_wall_v3_B = tic;
futs_v3 = cell(N_WORKERS, 1);
for w = 1:N_WORKERS
    futs_v3{w} = parfeval(pool, @(C) fast_coverage_simulator_function_v3(C, false, false, false, []), 0, Cfg);
end
for w = 1:N_WORKERS, wait(futs_v3{w}); end
t_v3_B = toc(t_wall_v3_B);
fprintf('  new v3: %.2f s wall\n', t_v3_B);

speedup_B_v2 = t_old_B / t_new_B;
speedup_B_v3 = t_old_B / t_v3_B;
fprintf('  >> v2 speedup: %.1fx  |  v3 speedup: %.1fx\n\n', speedup_B_v2, speedup_B_v3);

% =========================================================================
% SCENARIO C: Memory pressure sweep
%   Vary the number of simultaneous workers from 1 to N_WORKERS and record
%   wall time.  Shows where the L3 cache thrash knee is for each function.
%   sat_pos_ecef = 3 × num_sats × nT × 8 bytes
%   i9-13900K L3 = 36 MB, so thrash starts at 36 MB / sat_ecef_bytes workers.
% =========================================================================
sat_ecef_MB_d = 3 * Cfg.Total_sats * n_steps * 8 / 1e6;   % double
sat_ecef_MB_s = 3 * Cfg.Total_sats * n_steps * 4 / 1e6;   % single
thrash_at_d   = floor(36 / sat_ecef_MB_d);
thrash_at_s   = floor(36 / sat_ecef_MB_s);
fprintf('--- Scenario C: memory-pressure sweep ---\n');
fprintf('  sat_pos_ecef double = %.1f MB  →  L3 full at ~%d workers\n', sat_ecef_MB_d, thrash_at_d);
fprintf('  sat_pos_ecef single = %.1f MB  →  L3 full at ~%d workers\n\n', sat_ecef_MB_s, thrash_at_s);

% Also report P-core L2 fit (2 MB per P-core)
fprintf('  P-core L2 (2 MB):  double fits? %s   single fits? %s\n\n', ...
    string(sat_ecef_MB_d <= 2), string(sat_ecef_MB_s <= 2));

worker_counts = unique([1, 2, 4, 8, 16, N_WORKERS]);
worker_counts = worker_counts(worker_counts <= nw);

t_old_C = zeros(size(worker_counts));
t_new_C = zeros(size(worker_counts));
t_v3_C  = zeros(size(worker_counts));

for ci = 1:numel(worker_counts)
    npar = worker_counts(ci);

    futs = cell(npar, 1);
    tw = tic;
    for w = 1:npar
        futs{w} = parfeval(pool, @(C) coverage_simulator_function(C, false, false, []), 0, Cfg);
    end
    for w = 1:npar, wait(futs{w}); end
    t_old_C(ci) = toc(tw);

    futs = cell(npar, 1);
    tw = tic;
    for w = 1:npar
        futs{w} = parfeval(pool, @(C) fast_coverage_simulator_function_v2(C, false, false, false, []), 0, Cfg);
    end
    for w = 1:npar, wait(futs{w}); end
    t_new_C(ci) = toc(tw);

    futs = cell(npar, 1);
    tw = tic;
    for w = 1:npar
        futs{w} = parfeval(pool, @(C) fast_coverage_simulator_function_v3(C, false, false, false, []), 0, Cfg);
    end
    for w = 1:npar, wait(futs{w}); end
    t_v3_C(ci) = toc(tw);

    fprintf('  %2d workers:  old %.2f s  |  v2 %.2f s (%.1fx)  |  v3 %.2f s (%.1fx)\n', ...
        npar, t_old_C(ci), t_new_C(ci), t_old_C(ci)/t_new_C(ci), t_v3_C(ci), t_old_C(ci)/t_v3_C(ci));
end

% =========================================================================
% SUMMARY TABLE
% =========================================================================
fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('  %-28s  %8s  %8s  %8s  %7s  %7s\n', 'Scenario', 'Old (s)', 'v2 (s)', 'v3 (s)', 'v2 spd', 'v3 spd');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('  %-28s  %8.2f  %8.2f  %8.2f  %6.1fx  %6.1fx\n', 'A: single worker', t_old_A, t_new_A, t_v3_A, speedup_A_v2, speedup_A_v3);
fprintf('  %-28s  %8.2f  %8.2f  %8.2f  %6.1fx  %6.1fx\n', sprintf('B: %d workers (wall)', N_WORKERS), t_old_B, t_new_B, t_v3_B, speedup_B_v2, speedup_B_v3);
fprintf('%s\n', repmat('=', 1, 78));
fprintf('\n  L3 thrash:  double = ~%d workers  |  single = ~%d workers\n', thrash_at_d, thrash_at_s);
fprintf('  P-core L2 (2MB): double fits? %s  |  single fits? %s\n\n', ...
    string(sat_ecef_MB_d <= 2), string(sat_ecef_MB_s <= 2));

% =========================================================================
% FIGURE: wall-clock time vs number of simultaneous workers
% =========================================================================
c_old = [0.850, 0.325, 0.098];
c_new = [0.466, 0.674, 0.188];
c_v3  = [0.494, 0.184, 0.557];   % purple for v3 (single)

f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 720, 460]);
ax = axes(f1, 'Color', 'w');
hold(ax, 'on');

plot(ax, worker_counts, t_old_C, '-o', 'LineWidth', 2.5, 'Color', c_old, ...
    'MarkerFaceColor', c_old, 'MarkerSize', 8, 'DisplayName', 'original (double struct)');
plot(ax, worker_counts, t_new_C, '-s', 'LineWidth', 2.5, 'Color', c_new, ...
    'MarkerFaceColor', c_new, 'MarkerSize', 8, 'DisplayName', 'v2: matrix, double');
plot(ax, worker_counts, t_v3_C,  '-^', 'LineWidth', 2.5, 'Color', c_v3, ...
    'MarkerFaceColor', c_v3,  'MarkerSize', 8, 'DisplayName', 'v3: matrix, single (L2 fit)');

% Mark L3 thrash knees
xline(ax, thrash_at_d, '--', 'Color', c_new*0.7, 'LineWidth', 1.2, ...
    'Label', sprintf('double L3 full (%d)', thrash_at_d), ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(ax, thrash_at_s, '--', 'Color', c_v3*0.7,  'LineWidth', 1.2, ...
    'Label', sprintf('single L3 full (%d)', thrash_at_s), ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

xlabel(ax, 'Simultaneous workers',     'FontWeight', 'bold', 'FontSize', 13);
ylabel(ax, 'Wall-clock time (s)',       'FontWeight', 'bold', 'FontSize', 13);
legend(ax, 'Location', 'northwest',    'FontSize', 11);
grid(ax, 'on');  ax.Box = 'off';
title(ax, sprintf('Stage-3 equivalent: %d UEs \\times %d steps \\times %d sats', ...
    n_ues_stage3, n_steps, Cfg.Total_sats), ...
    'FontSize', 11, 'FontWeight', 'bold');

% Annotate speedup of v3 at 32 workers
x_end = max(worker_counts);
gap   = max([t_old_C(end), t_new_C(end), t_v3_C(end)]) * 0.04;
text(ax, x_end, t_new_C(end) + gap, sprintf('v2: %.0fx', t_old_C(end)/t_new_C(end)), ...
    'HorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold', 'Color', c_new);
text(ax, x_end, t_v3_C(end)  - gap*2, sprintf('v3: %.0fx', t_old_C(end)/t_v3_C(end)), ...
    'HorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold', 'Color', c_v3);

exportgraphics(f1, fullfile(out_dir, 'simulator_speed_comparison.png'), 'Resolution', 300);
close(f1);
fprintf('Saved: simulator_speed_comparison.png\n');

% =========================================================================
% FIGURE 2: time breakdown across all 3 functions (stacked bar)
%   v1: struct overhead (t_old_A - t_v2_A) + geometry + precision conversion (0)
%   v2: geometry only (double)
%   v3: geometry only (single) — expect noticeably lower due to L2 cache hits
%   Bars are single-worker times (Scenario A) to isolate overhead from contention.
% =========================================================================
t_struct_overhead = t_old_A - t_new_A;   % v1 - v2 = pure struct CoW cost
t_precision_gain  = t_new_A - t_v3_A;   % v2 - v3 = single-precision L2 benefit

% stacked: [struct_overhead, precision_overhead, geometry] for each function
bar_matrix = [ t_struct_overhead,  t_precision_gain,  t_v3_A;  % v1 total
               0,                  t_precision_gain,  t_v3_A;  % v2 total
               0,                  0,                 t_v3_A]; % v3 total

f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 440, 380]);
ax2 = axes(f2, 'Color', 'w');

bh = bar(ax2, bar_matrix, 'stacked');
bh(1).FaceColor = [0.850, 0.325, 0.098];   bh(1).DisplayName = 'Struct CoW overhead';
bh(2).FaceColor = [0.929, 0.694, 0.125];   bh(2).DisplayName = 'double → L3 spill (avoidable)';
bh(3).FaceColor = [0.18,  0.38,  0.60];    bh(3).DisplayName = 'Geometry (single, L2-resident)';

totals = sum(bar_matrix, 2)';
for k = 1:3
    text(ax2, k, totals(k)*1.03, sprintf('%.1f s', totals(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end

ax2.XTickLabel = {'v1 original', 'v2 (matrix)', 'v3 (single)'};
ylabel(ax2, 'Single-worker time (s)', 'FontWeight', 'bold', 'FontSize', 13);
legend(ax2, 'Location', 'northeast', 'FontSize', 10);
grid(ax2, 'on');  ax2.Box = 'off';  ax2.XGrid = 'off';
title(ax2, sprintf('Cost breakdown  |  %d UEs \\times %d steps  |  %d sats', ...
    n_ues_stage3, n_steps, Cfg.Total_sats), 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f2, fullfile(out_dir, 'simulator_struct_overhead.png'), 'Resolution', 300);
close(f2);
fprintf('Saved: simulator_struct_overhead.png\n');
fprintf('\nAll done.\n');
