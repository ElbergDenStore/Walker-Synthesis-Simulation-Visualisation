% PLOT_GRIDSEARCH  Regenerate gridsearch plots from saved plot_data.mat files.
%
% Set TARGET below, then run:
%   ''                  - most recent Master_Sweep_*/gridsearch_runs/
%   'path/to/sweep'     - all runs inside a Master_Sweep folder
%   'path/to/run'       - one specific gridsearch run folder
%   {'path1','path2'}   - list of run folders
target = '/home/aau/master_matlab_sims/simulation_output/gridsearch_runs/8000_18_IRIS2GlobalMEO';
% ├── 1200_264_IRIS2GLOBALDELTA
% │   ├── deep_profile_report.txt
% │   └── plot_data.mat
% ├── 1200_264_IRIS2GLOBALSTAR
% │   ├── deep_profile_report.txt
% │   └── plot_data.mat
% ├── 1200_264_IRIS2REGIONALDELTA
% │   ├── deep_profile_report.txt
% │   └── plot_data.mat
% ├── 1200_264_IRIS2REGIONALDELTAOLD
% /home/aau/master_matlab_sims/simulation_output/gridsearch_runs/8000_18_IRIS2GlobalMEO

% Add the project to the MATLAB path (robust to the script's folder depth).
workspace_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(workspace_root, 'functions', 'path_setup.m')), workspace_root = fileparts(workspace_root); end
addpath(fullfile(workspace_root, 'functions'));
path_setup();

sim_dir   = fullfile(workspace_root, 'simulation_output');
runs_root = fullfile(sim_dir, 'gridsearch_runs');

if isempty(target)
    % Prefer the most recent Master_Sweep folder's nested gridsearch_runs;
    % fall back to the legacy top-level gridsearch_runs directory.
    sweep_hits = dir(fullfile(sim_dir, 'Master_Sweep_*'));
    sweep_hits = sweep_hits([sweep_hits.isdir]);
    if ~isempty(sweep_hits)
        [~, best] = max([sweep_hits.datenum]);
        search_root = fullfile(sim_dir, sweep_hits(best).name, 'gridsearch_runs');
        fprintf('Searching within most recent sweep: %s\n', search_root);
    else
        search_root = runs_root;
    end
    hits = dir(fullfile(search_root, '*', 'plot_data.mat'));
    if isempty(hits)
        error('No plot_data.mat files found under:\n  %s', search_root);
    end
    folders = unique({hits.folder});
    fprintf('Found %d run(s) with plot data.\n', numel(folders));
elseif ischar(target) || isstring(target)
    target = char(target);
    if target(1) ~= '/' && ~(numel(target) > 1 && target(2) == ':')
        target = fullfile(workspace_root, target);
    end
    % If target is a Master_Sweep folder, expand to all nested runs.
    nested_hits = dir(fullfile(target, 'gridsearch_runs', '*', 'plot_data.mat'));
    if ~isempty(nested_hits)
        folders = unique({nested_hits.folder});
        fprintf('Master_Sweep folder: found %d run(s) with plot data.\n', numel(folders));
    else
        folders = {target};
    end
elseif iscell(target)
    folders = target;
else
    error('Input must be a folder path (string) or cell array of paths.');
end

fprintf('Regenerating gridsearch plots for %d run(s)...\n', numel(folders));
for i = 1:numel(folders)
    folder   = abs_path(folders{i}, workspace_root);
    mat_file = fullfile(folder, 'plot_data.mat');
    if ~isfile(mat_file)
        warning('No plot_data.mat in:\n  %s\nSkipping.', folder);
        continue;
    end
    fprintf('  [%d/%d] %s\n', i, numel(folders), folder);
    s = load(mat_file, 'plot_data');
    generate_plots(s.plot_data, folder);
end
fprintf('Done.\n');

% -------------------------------------------------------------------------
function generate_plots(pd, out_dir)

search_grid           = pd.search_grid;
status_flags          = pd.status_flags;
evaluated_coverage    = pd.evaluated_coverage;
detailed_coverage     = pd.detailed_coverage;
orbit_height_km       = pd.orbit_height_km;
worst_accepted_sats   = pd.worst_accepted_sats;
all_candidates        = pd.all_candidates; %#ok<NASGU>
target_num_candidates = pd.target_num_candidates; %#ok<NASGU>

%% Styling
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);
color_inv  = [0.8 0.8 0.8];
color_cand = [0.2 0.6 0.8];
color_best = [1.0 0.8 0.0];
sz_inv  = 35;  sz_cand = 50;  sz_best = 120;
export_dpi = 300;
leg_loc    = 'northeast';
fig_pos    = [100, 100, 500, 400];

%% Derived data
was_evaluated = ~isnan(evaluated_coverage);
eval_grid     = search_grid(was_evaluated, :);
eval_status   = status_flags(was_evaluated);
isInvalid = eval_status == 0;
isCand    = eval_status == 1;
isBest    = eval_status == 2;

history_Loss = eval_grid.Total_sats;
phasing_deg  = (eval_grid.Phasing ./ eval_grid.Num_planes) .* 360;
history_X    = table(eval_grid.Num_planes, eval_grid.Sats_per_plane, ...
    eval_grid.Inclination, phasing_deg, ...
    'VariableNames', {'Num_planes', 'Sats_per_plane', 'Inclination', 'Phasing_Degrees'});

% Shared colour limits (num sats) — consistent across all plots
cand_or_best = history_Loss(isCand | isBest);
if isempty(cand_or_best) || min(cand_or_best) == max(cand_or_best)
    clim_val = [min(history_Loss)-1, max(history_Loss)+1];
else
    clim_val = [min(cand_or_best), max(cand_or_best)];
end

%% Plot 1: Num Sats vs Inclination
f1 = figure('Visible', 'off', 'Name', 'Sats vs Inclination', 'Color', 'w', 'Position', fig_pos);
hold on;
[jx_inc, jy_loss] = apply_density_jitter(history_X.Inclination, history_Loss, false, true);
scatter(jx_inc(isInvalid), jy_loss(isInvalid), sz_inv,  color_inv,  'x');
scatter(jx_inc(isCand),    jy_loss(isCand),    sz_cand, history_Loss(isCand), 'filled', 'MarkerEdgeColor', 'k');
scatter(jx_inc(isBest),    jy_loss(isBest),    sz_best, history_Loss(isBest), 'diamond', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
if ~isnan(worst_accepted_sats)
    yl_val = worst_accepted_sats + 0.5;
    yline(yl_val, '--r', 'LineWidth', 1.5, 'Color', [0.8 0 0]);
    xl = xlim;
    text(xl(1), yl_val, ' Exhausted Search Space Below', ...
        'Color', [0.8 0 0], 'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'white', 'EdgeColor', 'none', 'Margin', 3);
end
xlabel('Inclination (deg)', 'FontWeight', 'bold'); ylabel('Satellite Count', 'FontWeight', 'bold');
% title("Candidate Inclinations @ " + num2str(orbit_height_km) + " km");
legend('Invalid', 'Candidate', 'Minimum', 'Location', leg_loc);
colormap(f1, 'parula'); set(gca, 'CLim', clim_val);
grid on; hold off;
exportgraphics(f1, fullfile(out_dir, 'Inclinations_NumSats.png'), 'Resolution', export_dpi);
close(f1);

%% Plot 1b: Coverage vs Total Sats (detailed runs only)
was_detailed = ~isnan(detailed_coverage);
if any(was_detailed)
    det_grid   = search_grid(was_detailed, :);
    det_cov    = detailed_coverage(was_detailed);
    det_status = status_flags(was_detailed);

    f_trade = figure('Visible', 'off', 'Name', 'Detailed Tradeoff', 'Color', 'w', 'Position', fig_pos);
    hold on;
    [jx_sats, jy_cov] = apply_density_jitter(det_grid.Total_sats, det_cov, true, false);
    scatter(jx_sats(det_status==0), jy_cov(det_status==0), sz_inv,  color_inv,  'x', 'LineWidth', 1.0);
    scatter(jx_sats(det_status==1), jy_cov(det_status==1), sz_cand, det_grid.Total_sats(det_status==1), 'filled', 'MarkerEdgeColor', 'k');
    scatter(jx_sats(det_status==2), jy_cov(det_status==2), sz_best, det_grid.Total_sats(det_status==2), 'diamond', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    xlabel('Satellite Count', 'FontWeight', 'bold'); ylabel('Worst Coverage %', 'FontWeight', 'bold');
    % title("Coverage percentage @ " + num2str(orbit_height_km) + " km");
    legend('Invalid', 'Candidate', 'Minimum', 'Location', leg_loc);
    colormap(f_trade, 'parula'); set(gca, 'CLim', clim_val);
    grid on; hold off;
    exportgraphics(f_trade, fullfile(out_dir, 'Detailed_Tradeoff_Coverage_vs_Sats.png'), 'Resolution', export_dpi);
    close(f_trade);
end

%% Plot 2: Architecture map (planes vs sats per plane)
planes  = history_X.Num_planes;
sats_pp = history_X.Sats_per_plane;

f4 = figure('Visible', 'off', 'Name', 'Architecture Map', 'Color', 'w', 'Position', fig_pos);
hold on;
[jx_planes, jy_sats_pp] = apply_density_jitter(planes, sats_pp, true, true);
scatter(jx_planes(isInvalid), jy_sats_pp(isInvalid), sz_inv,  color_inv,  'x');
scatter(jx_planes(isCand),    jy_sats_pp(isCand),    sz_cand, history_Loss(isCand), 'filled', 'MarkerEdgeColor', 'k');
scatter(jx_planes(isBest),    jy_sats_pp(isBest),    sz_best, history_Loss(isBest), 'diamond', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
xlabel('Orbital Plane Count', 'FontWeight', 'bold'); ylabel('Satellites per Plane', 'FontWeight', 'bold');
% title("Evaluated Constellations @ " + num2str(orbit_height_km) + " km");
legend('Invalid', 'Candidate', 'Minimum', 'Location', leg_loc);
colormap(f4, 'parula'); set(gca, 'CLim', clim_val);
grid on; hold off;
exportgraphics(f4, fullfile(out_dir, 'NumPlanes_SatsPerPlane.png'), 'Resolution', export_dpi);
close(f4);

%% Plot 3: Phasing vs Number of Planes
f_phase = figure('Visible', 'off', 'Name', 'Phasing vs Planes', 'Color', 'w', 'Position', fig_pos);
hold on;
[jx_phase_p, jy_phase_deg] = apply_density_jitter(history_X.Num_planes, history_X.Phasing_Degrees, true, false);
scatter(jx_phase_p(isInvalid), jy_phase_deg(isInvalid), sz_inv,  color_inv,  'x');
scatter(jx_phase_p(isCand),    jy_phase_deg(isCand),    sz_cand, history_Loss(isCand), 'filled', 'MarkerEdgeColor', 'k');
scatter(jx_phase_p(isBest),    jy_phase_deg(isBest),    sz_best, history_Loss(isBest), 'diamond', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
xlabel('Orbital Plane Count', 'FontWeight', 'bold'); ylabel('Normalised Phasing (deg)', 'FontWeight', 'bold');
% title("Candidate Phasing @ " + num2str(orbit_height_km) + " km");
legend('Invalid', 'Candidate', 'Minimum', 'Location', leg_loc);
ylim([0 360]); yticks(0:45:360);
colormap(f_phase, 'parula'); set(gca, 'CLim', clim_val);
cb = colorbar; cb.Label.String = 'Satellite Count';
grid on; hold off;
exportgraphics(f_phase, fullfile(out_dir, 'NumPlanes_Phasing.png'), 'Resolution', export_dpi);
close(f_phase);

fprintf('    Saved 4 plots to: %s\n', out_dir);
end

% -------------------------------------------------------------------------
function p = abs_path(p, root)
% Return absolute path, resolving relative paths against root.
if isempty(p), return; end
if p(1) ~= '/' && ~(numel(p) > 1 && p(2) == ':')
    p = fullfile(root, p);
end
end

function [jx, jy] = apply_density_jitter(x, y, jitter_x_flag, jitter_y_flag)
jx = x; jy = y;
[~, ~, ic] = unique([x, y], 'rows');
counts = accumarray(ic, 1);
for i = 1:max(ic)
    n = counts(i);
    if n > 1
        idx    = find(ic == i);
        spread = min(0.20, 0.02 * sqrt(n - 1));
        if jitter_x_flag, jx(idx) = x(idx) + (rand(length(idx), 1) - 0.5) * 2 * spread; end
        if jitter_y_flag, jy(idx) = y(idx) + (rand(length(idx), 1) - 0.5) * 2 * spread; end
    end
end
end
