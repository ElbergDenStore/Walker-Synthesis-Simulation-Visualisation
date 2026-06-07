function ok = smoke_test(verbose)
% SMOKE_TEST  Fast safety net for the repo cleanup (renames / file moves).
%
%   ok = smoke_test()        % run all checks, print a summary
%   ok = smoke_test(true)    % also print per-file checkcode detail
%
% Runs three layers, cheapest first:
%
%   1. STATIC  - checkcode() over every tracked *.m file (excludes old/).
%                Catches syntax/parse errors introduced by edits. NOTE:
%                MATLAB's analyzer cannot flag a call to a *renamed* function
%                as undefined, so this layer will NOT catch dangling callers
%                on its own - layers 2/3 and the rename-map grep do that.
%
%   2. RUNTIME - builds a tiny config via default_config() and runs the full
%                Constellation_simulator() pipeline (generators -> beams ->
%                geometry -> coverage). This DOES throw if a function in the
%                core call graph was renamed/moved and a caller wasn't updated.
%
%   3. DEPS    - matlab.codetools.requiredFilesAndProducts() on each entry
%                point. Resolves the dependency graph; a parse/availability
%                problem surfaces here as an error.
%
% Returns true only if every layer passes. Designed to run headless
% (no plotting, no parallel pool, no Satellite Toolbox propagator).
%
% Recommended use during cleanup: run before a batch of renames to get a
% green baseline, then again after, and diff the result.

    if nargin < 1, verbose = false; end

    % ---- Path bootstrap (independent of caller cwd) ---------------------
    % Add the project to the MATLAB path (robust to the script's folder depth).
    repo_root = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
    addpath(fullfile(repo_root, 'functions'));
    path_setup();

    fprintf('\n===== SMOKE TEST =====\n');
    fprintf('Repo: %s\n', repo_root);

    results = struct('name', {}, 'pass', {}, 'detail', {});

    % =====================================================================
    % LAYER 1 - static checkcode sweep
    % =====================================================================
    fprintf('\n[1/3] Static checkcode sweep...\n');
    files = dir(fullfile(repo_root, '**', '*.m'));
    n_err_files = 0;
    n_scanned   = 0;
    for k = 1:numel(files)
        fpath = fullfile(files(k).folder, files(k).name);
        rel   = strrep(fpath, [repo_root filesep], '');

        % Skip the archive and this harness itself.
        if startsWith(rel, ['old' filesep])
            continue;
        end
        n_scanned = n_scanned + 1;

        msgs = checkcode(fpath, '-id');
        % Keep only genuine parse/syntax errors, not style warnings.
        is_err = arrayfun(@(m) is_syntax_error(m), msgs);
        errs   = msgs(is_err);

        if ~isempty(errs)
            n_err_files = n_err_files + 1;
            fprintf('  ERROR  %s\n', rel);
            for e = errs(:)'
                fprintf('         line %d: %s\n', e.line, e.message);
            end
        elseif verbose
            fprintf('  ok     %s\n', rel);
        end
    end
    pass1 = (n_err_files == 0);
    results(end+1) = mk('checkcode', pass1, ...
        sprintf('%d files scanned, %d with syntax errors', n_scanned, n_err_files));

    % =====================================================================
    % LAYER 2 - tiny end-to-end run through the core pipeline
    % =====================================================================
    fprintf('\n[2/3] Runtime smoke (Constellation_simulator, tiny config)...\n');
    pass2  = true;
    detail2 = '';
    try
        % Smallest meaningful config: ~6 UEs, 1 hour, fast-math propagator.
        Cfg = default_config(1000, "walkerdelta", "small", "short");
        metrics = Constellation_simulator(Cfg, false, false, false);  % no parallel/link/toolbox

        assert(isstruct(metrics), 'metrics is not a struct');
        assert(isfield(metrics, 'worst_coverage_percent'), 'missing worst_coverage_percent');
        wc = metrics.worst_coverage_percent;
        assert(isscalar(wc) && isfinite(wc) && wc >= 0 && wc <= 100, ...
            sprintf('worst_coverage_percent out of range: %g', wc));
        detail2 = sprintf('worst_coverage_percent = %.2f%%', wc);
        fprintf('  ok     %s\n', detail2);
    catch ME
        pass2  = false;
        detail2 = sprintf('%s: %s', ME.identifier, ME.message);
        fprintf('  ERROR  %s\n', detail2);
    end
    results(end+1) = mk('runtime', pass2, detail2);

    % =====================================================================
    % LAYER 3 - dependency resolution on entry points
    % =====================================================================
    fprintf('\n[3/3] Dependency resolution on entry points...\n');
    entry_points = {
        'Constellation_simulator.m'
        'run_constellation_simulation_example.m'
        'numerical_walker_synthesis.m'
        'analytical_walker_star_synthesis.m'
        'plot_simulation.m'
        fullfile('functions', 'gridsearch.m')
        'default_config.m'
    };
    pass3 = true;
    for k = 1:numel(entry_points)
        ep = fullfile(repo_root, entry_points{k});
        if ~isfile(ep)
            fprintf('  skip   %s (not found)\n', entry_points{k});
            continue;
        end
        try
            req = matlab.codetools.requiredFilesAndProducts(ep);
            fprintf('  ok     %-34s (%d files)\n', entry_points{k}, numel(req));
        catch ME
            pass3 = false;
            fprintf('  ERROR  %-34s %s\n', entry_points{k}, ME.message);
        end
    end
    results(end+1) = mk('dependencies', pass3, '');

    % =====================================================================
    % Summary
    % =====================================================================
    ok = all([results.pass]);
    fprintf('\n----- SUMMARY -----\n');
    for r = results
        fprintf('  %-14s : %s%s\n', r.name, tern(r.pass, 'PASS', 'FAIL'), ...
            tern(isempty(r.detail), '', [' (' r.detail ')']));
    end
    fprintf('  %-14s : %s\n', 'OVERALL', tern(ok, 'PASS', 'FAIL'));
    fprintf('===================\n\n');
end

% ------------------------------------------------------------------------
function tf = is_syntax_error(m)
    % checkcode message ids for hard parse/syntax errors generally start with
    % these prefixes; treat anything mentioning a parse failure as an error.
    tf = contains(lower(m.message), 'parse error') || ...
         contains(lower(m.message), 'invalid syntax') || ...
         startsWith(m.id, 'MATLAB:lang') || ...
         strcmp(m.id, 'pZZ');   % generic "unable to parse" id
end

function s = mk(name, pass, detail)
    s = struct('name', name, 'pass', logical(pass), 'detail', detail);
end

function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end
