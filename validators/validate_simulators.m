function ok = validate_simulators(verbose)
% VALIDATE_SIMULATORS  Prove constellation_simulator is correct across modes.
%
%   ok = validate_simulators()        % run all checks, print summary
%   ok = validate_simulators(true)    % print extra per-test detail
%
% constellation_simulator has two propagation backends, selected by the
% use_toolbox flag:
%
%   use_toolbox = true   Satellite Toolbox states(), driven through a
%                        PERSISTENT cached satelliteScenario (see
%                        reuse_or_create_scenario in constellation_simulator).
%                        The scenario container is reused between calls; only
%                        satellites/ground stations are rebuilt. This reuse is
%                        the "caching".
%   use_toolbox = false  Pure-math fast_walker_ecef (Walker Star + Delta).
%                        No toolbox handle objects, no persistent state.
%
% The deterministic geometry output is metrics.Num_visible (numUEs x nT,
% the count of satellites above Min_elevation_UE for each UE at each step).
% Every test below compares that matrix.
%
% TEST 1 - CACHE INTEGRITY (does caching ruin anything?)
%   Run config A cold (cleared cache) -> A1.
%   Run a DIFFERENT config B          -> pollutes the persistent scenario.
%   Run config A again                -> A2.
%   If the cache leaks any state between calls, A2 differs from A1.
%   PASS requires A1 and A2 to be BIT-IDENTICAL. A back-to-back repeat of A
%   is also checked for exact reproducibility.
%
% TEST 2 - FAST vs TOOLBOX EQUIVALENCE
%   Run the same config A through both backends and compare. They use
%   independent code paths (different propagators), so agreement is checked
%   within a small tolerance rather than bit-for-bit.
%
% Returns true only if every test passes. Runs headless: tiny UE grid,
% short duration, no parallel pool, no link budget.

    if nargin < 1, verbose = false; end

    % ---- Path bootstrap (independent of caller cwd) ---------------------
    repo_root = fileparts(fileparts(mfilename('fullpath')));  % validators/ -> root
    addpath(fullfile(repo_root, 'functions'));
    addpath(repo_root);
    path_setup();

    % Tolerances for TEST 2 (different propagators, not bit-identical).
    COV_TOL_PP   = 1.0;   % max allowed |worst_coverage_percent| difference [pp]
    MEAN_VIS_TOL = 0.05;  % max allowed mean |visible-count| difference [sats]

    fprintf('\n===== VALIDATE SIMULATORS =====\n');

    % Two genuinely different constellations so the cache is really stressed:
    % different altitude, type, and geometry.
    cfgA = get_cfg(1000, "walkerdelta", "small", "short");
    cfgB = get_cfg( 600, "walkerstar",  "small", "short");
    fprintf('Config A: walkerdelta 1000 km, %d sats (%dx%d)\n', ...
        cfgA.Total_sats, cfgA.Num_planes, cfgA.Sats_per_plane);
    fprintf('Config B: walkerstar   600 km, %d sats (%dx%d)\n', ...
        cfgB.Total_sats, cfgB.Num_planes, cfgB.Sats_per_plane);

    results = struct('name', {}, 'pass', {}, 'detail', {});

    % =====================================================================
    % TEST 1 - cache integrity
    % =====================================================================
    fprintf('\n[1] Cache integrity (toolbox backend)...\n');

    clear constellation_simulator;   % reset the persistent cached scenario

    A1  = run_q(cfgA, true);         % A cold
    A1b = run_q(cfgA, true);         % A again, back-to-back (warm, same cfg)
    run_q(cfgB, true);               % pollute cache with a different B
    A2  = run_q(cfgA, true);         % A after a different constellation ran

    repeat_ok  = isequal(A1, A1b);
    pollute_ok = isequal(A1, A2);
    pass1 = repeat_ok && pollute_ok;

    if verbose || ~repeat_ok
        fprintf('   back-to-back repeat identical    : %s\n', yn(repeat_ok));
    end
    if verbose || ~pollute_ok
        fprintf('   identical after B polluted cache : %s\n', yn(pollute_ok));
    end
    if ~pollute_ok
        d = abs(A1 - A2);
        fprintf('   (cache leak) cells differing: %d, max diff: %g\n', ...
            nnz(d), max(d(:)));
    end
    results(end+1) = mk('cache_integrity', pass1, ...
        sprintf('repeat=%s, post-pollute=%s', yn(repeat_ok), yn(pollute_ok)));

    % =====================================================================
    % TEST 2 - fast vs toolbox equivalence
    % =====================================================================
    fprintf('\n[2] Fast vs toolbox equivalence (config A)...\n');

    M_tb   = run_full(cfgA, true);    % toolbox (cached scenario)
    M_fast = run_full(cfgA, false);   % pure-math

    % The two backends build their own time grids and can differ by one
    % endpoint sample (toolbox may include the StopTime instant). Both start
    % at StartTime with the same SampleTime, so the first nT_common columns
    % are the same timestamps; compare on that aligned window.
    V_tb   = double(M_tb.Num_visible);
    V_fast = double(M_fast.Num_visible);
    nT_tb   = size(V_tb, 2);
    nT_fast = size(V_fast, 2);
    nT_common = min(nT_tb, nT_fast);
    if nT_tb ~= nT_fast
        fprintf('   note: time steps differ (toolbox=%d, fast=%d); comparing first %d\n', ...
            nT_tb, nT_fast, nT_common);
    end
    V_tb   = V_tb(:,   1:nT_common);
    V_fast = V_fast(:, 1:nT_common);

    d         = abs(V_tb - V_fast);
    max_diff  = max(d(:));
    mean_diff = mean(d(:));
    pct_cells = 100 * nnz(d) / numel(d);
    cov_diff  = abs(M_tb.worst_coverage_percent - M_fast.worst_coverage_percent);

    pass2 = (cov_diff <= COV_TOL_PP) && (mean_diff <= MEAN_VIS_TOL);

    fprintf('   worst_coverage_percent : toolbox=%.3f  fast=%.3f  |diff|=%.3f pp (tol %.2f)\n', ...
        M_tb.worst_coverage_percent, M_fast.worst_coverage_percent, cov_diff, COV_TOL_PP);
    fprintf('   visible-count matrix   : mean|diff|=%.4f (tol %.2f)  max|diff|=%g  cells differing=%.2f%%\n', ...
        mean_diff, MEAN_VIS_TOL, max_diff, pct_cells);
    results(end+1) = mk('fast_vs_toolbox', pass2, ...
        sprintf('cov|diff|=%.3f pp, mean|diff|=%.4f sats', cov_diff, mean_diff));

    % =====================================================================
    % Summary
    % =====================================================================
    ok = all([results.pass]);
    fprintf('\n----- SUMMARY -----\n');
    for r = results
        fprintf('  %-16s : %s (%s)\n', r.name, yn(r.pass), r.detail);
    end
    fprintf('  %-16s : %s\n', 'OVERALL', yn(ok));
    fprintf('===============================\n\n');
end

% ------------------------------------------------------------------------
function V = run_q(Cfg, use_toolbox)
% Quiet run, returns only the deterministic visibility matrix (double).
    m = run_full(Cfg, use_toolbox);
    V = double(m.Num_visible);
end

function m = run_full(Cfg, use_toolbox)
% Run the simulator with stdout suppressed (keeps the report clean).
    evalc('m = constellation_simulator(Cfg, false, false, use_toolbox);');
end

function s = mk(name, pass, detail)
    s = struct('name', name, 'pass', logical(pass), 'detail', detail);
end

function s = yn(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end
