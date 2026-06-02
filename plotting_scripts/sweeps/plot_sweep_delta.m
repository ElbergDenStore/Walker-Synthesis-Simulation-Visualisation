function plot_sweep_delta(sweep_folder)
% PLOT_SWEEP_DELTA  Plot and summarize a Walker Delta gridsearch sweep.
%   Compares the analytical Walker-Star formula (SOC baseline) against
%   the numerical Walker-Delta gridsearch result at each altitude.
%
%   Outputs (saved to sweep folder):
%     Star_vs_Delta_Comparison.png   – scatter comparison plot
%     walker_comparison_table.tex    – IEEE LaTeX table
%     multistage_table.tex           – multi-stage filtering stats
%
% Usage:
%   plot_sweep_delta()                  – most recent Master_Sweep_* (non-WalkerStar) folder
%   plot_sweep_delta('simulation_output/Master_Sweep_20260522_172327')   – specific sweep folder

script_dir     = fileparts(mfilename('fullpath'));
workspace_root = fileparts(fileparts(script_dir));
sim_dir        = fullfile(workspace_root, 'simulation_output');

%% ---- Locate sweep folder -----------------------------------------------
if nargin == 0
    hits = dir(fullfile(sim_dir, 'Master_Sweep_*'));
    hits = hits([hits.isdir]);
    % Exclude Walker Star sweeps so we pick up delta-only runs
    hits = hits(~contains({hits.name}, 'WalkerStar'));
    if isempty(hits)
        error('No Master_Sweep_* (delta) folders found in:\n  %s', sim_dir);
    end
    [~, idx]     = max([hits.datenum]);
    sweep_folder = fullfile(sim_dir, hits(idx).name);
    fprintf('Using most recent Delta sweep: %s\n', sweep_folder);
else
    sweep_folder = char(sweep_folder);
    if sweep_folder(1) ~= '/' && ~(numel(sweep_folder) > 1 && sweep_folder(2) == ':')
        sweep_folder = fullfile(workspace_root, sweep_folder);
    end
end

mat_file = fullfile(sweep_folder, 'Master_Altitude_Sweep_Results.mat');
if ~isfile(mat_file)
    error('Master_Altitude_Sweep_Results.mat not found in:\n  %s', sweep_folder);
end

loaded = load(mat_file);
if ~isfield(loaded, 'best_delta_sats')
    error(['This folder contains a Walker Star sweep, not a Walker Delta sweep.\n' ...
           'Use plot_sweep_star() instead.']);
end

heights_km      = loaded.heights_km;
star_sats       = loaded.star_sats;       % analytical Walker Star (SOC baseline)
best_delta_sats = loaded.best_delta_sats; % numerical Walker Delta

% Sort ascending for consistent indexing
[heights_km, sort_idx] = sort(heights_km, 'ascend');
star_sats       = star_sats(sort_idx);
best_delta_sats = best_delta_sats(sort_idx, :);

star_T  = [star_sats.Total_sats]';
delta_T = best_delta_sats.Total_sats;

%% ---- Plot ---------------------------------------------------------------
if isfield(loaded, 'Master_config')
    cfg             = loaded.Master_config;
    minimum_lat_deg = min(cfg.Lat_range_deg);
    maximum_lat_deg = max(cfg.Lat_range_deg);
    title_str = sprintf( ...
        'Optimal Walker Constellations | \\lambda: %.1f^{\\circ} - %.1f^{\\circ} | \\epsilon_{min}: %.1f^{\\circ}', ...
        minimum_lat_deg, maximum_lat_deg, cfg.Min_elevation_UE);
else
    title_str = 'Walker Star (Analytical) vs Walker Delta (Numerical)';
end

set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

f1 = figure('Visible', 'off', 'Name', 'Constellation Comparison', 'Color', 'w', 'Position', [100 100 700 450]);
hold on;
scatter(heights_km, star_T,  36, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'Walker Star (Analytical)');
scatter(heights_km, delta_T, 36, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Walker Delta (Numerical)');
xlabel('Equatorial Orbital Altitude (km)', 'FontWeight', 'bold');
ylabel('Satellite Count', 'FontWeight', 'bold');
title(title_str);
legend('Location', 'northeast');
grid on; hold off;

plot_path = fullfile(sweep_folder, 'Star_vs_Delta_Comparison.png');
exportgraphics(f1, plot_path, 'Resolution', 300);
close(f1);
fprintf('Plot saved to: %s\n', plot_path);

%% ---- Parameter distribution plots -------------------------------------
generate_parameter_sweep_plots(heights_km, best_delta_sats, sweep_folder);

%% ---- Console summary ---------------------------------------------------
n          = numel(heights_km);
abs_saving = star_T - delta_T;
pct_saving = 100 * abs_saving ./ star_T;

fprintf('\n===== Walker-Star vs Walker-Delta Savings Summary =====\n');
fprintf('%-10s  %-8s  %-8s  %-10s  %-10s\n', ...
    'Alt (km)', 'Star-T', 'Delta-T', 'Abs Save', 'Pct Save (%)');
fprintf('%s\n', repmat('-', 1, 54));
for k = 1:n
    fprintf('%-10d  %-8d  %-8d  %-10d  %-10.1f\n', ...
        heights_km(k), star_T(k), delta_T(k), abs_saving(k), pct_saving(k));
end

fprintf('\nOverall range (all altitudes):\n');
fprintf('  Abs savings : %d – %d satellites\n', min(abs_saving), max(abs_saving));
fprintf('  Pct savings : %.1f%% – %.1f%%\n',    min(pct_saving), max(pct_saving));
if any(pct_saving < 0)
    neg_alts = heights_km(pct_saving < 0);
    fprintf('  NOTE: %d altitude(s) where Delta is worse than Star:\n', numel(neg_alts));
    for ki = 1:numel(neg_alts)
        ni = find(heights_km == neg_alts(ki), 1);
        fprintf('    %d km  (Star=%d, Delta=%d) — expected: Delta inclination capped at 80°\n', ...
            neg_alts(ki), star_T(ni), delta_T(ni));
    end
end

%% ---- Representative altitude detail -----------------------------------
rep_alts = [500, 700, 1000, 1200];
rep_idx  = zeros(1, numel(rep_alts));
for k = 1:numel(rep_alts)
    [~, rep_idx(k)] = min(abs(heights_km - rep_alts(k)));
end

fprintf('\n===== Representative Altitude Detail =====\n');
hdr = sprintf('%-6s | %-4s %-4s %-4s %-4s %-5s | %-4s %-4s %-4s %-4s %-5s | %-8s %-8s', ...
    'Alt', 'T*', 'P*', 'S*', 'F*', 'i*(°)', ...
    'T∆', 'P∆', 'S∆', 'F∆', 'i∆(°)', ...
    'ΔT', 'Save%');
fprintf('%s\n%s\n', hdr, repmat('-', 1, numel(hdr)));
for k = 1:numel(rep_idx)
    ri = rep_idx(k);
    ss = star_sats(ri);
    ds = best_delta_sats(ri, :);
    fprintf('%-6d | %-4d %-4d %-4d %-4d %-5.1f | %-4d %-4d %-4d %-4d %-5.2f | %-8d %-8.1f\n', ...
        heights_km(ri), ...
        ss.Total_sats, ss.Num_planes, ss.Sats_per_plane, ss.Phasing, ss.Inclination, ...
        ds.Total_sats, ds.Num_planes, ds.Sats_per_plane, ds.Phasing, ds.Inclination, ...
        ss.Total_sats - ds.Total_sats, ...
        100 * (ss.Total_sats - ds.Total_sats) / ss.Total_sats);
end

%% ---- LaTeX table (IEEE Access style) -----------------------------------
tex_lines = {};
tex_lines{end+1} = '% Auto-generated by plot_sweep_delta.m';
tex_lines{end+1} = '% Paste into your IEEE Access manuscript.';
tex_lines{end+1} = '%';
tex_lines{end+1} = '% Column key:';
tex_lines{end+1} = '%   T = total satellites, P = planes, S = sats/plane,';
tex_lines{end+1} = '%   F = phasing parameter, i = inclination (degrees)';
tex_lines{end+1} = '';
tex_lines{end+1} = '\begin{table}[!t]';
tex_lines{end+1} = '  \caption{Walker Constellation Comparison for Subset of Orbital Altitudes}';
tex_lines{end+1} = '  \label{tab:walker_comparison}';
tex_lines{end+1} = '  \centering';
tex_lines{end+1} = '  \renewcommand{\arraystretch}{1.2}';
tex_lines{end+1} = '  \begin{tabular}{c|cccc|ccccc|cc}';
tex_lines{end+1} = '    \hline\hline';
tex_lines{end+1} = '    & \multicolumn{4}{c|}{\textbf{Walker-Star (SOC)}}';
tex_lines{end+1} = '    & \multicolumn{5}{c|}{\textbf{Walker-Delta}}';
tex_lines{end+1} = '    & \multicolumn{2}{c}{\textbf{Difference}} \\';
tex_lines{end+1} = '    \textbf{Alt.} & $T$ & $P$ & $S$ & $i$ (°)';
tex_lines{end+1} = '    & $T$ & $P$ & $S$ & $F$ & $i$ (°)';
tex_lines{end+1} = '    & $\Delta T$ & \% \\';
tex_lines{end+1} = '    \textbf{(km)} & & & & & & & & & & & \\';
tex_lines{end+1} = '    \hline';

for k = 1:numel(rep_idx)
    ri  = rep_idx(k);
    ss  = star_sats(ri);
    ds  = best_delta_sats(ri, :);
    dT  = ss.Total_sats - ds.Total_sats;
    pct = 100 * dT / ss.Total_sats;
    tex_lines{end+1} = sprintf( ...   %#ok<AGROW>
        '    %d & %d & %d & %d & %.0f & %d & %d & %d & %d & %.2f & %d & %.1f\\%% \\\\', ...
        heights_km(ri), ...
        ss.Total_sats, ss.Num_planes, ss.Sats_per_plane, ss.Inclination, ...
        ds.Total_sats, ds.Num_planes, ds.Sats_per_plane, ds.Phasing, ds.Inclination, ...
        dT, pct);
end

tex_lines{end+1} = '    \hline\hline';
tex_lines{end+1} = '  \end{tabular}';
tex_lines{end+1} = '\end{table}';

%% ---- Abstract / conclusion sentence ------------------------------------
pct_min = min(pct_saving);
pct_max = max(pct_saving);
alt_min = min(heights_km);
alt_max = max(heights_km);
[r, ~]  = corrcoef(heights_km, pct_saving);
if r(1, 2) > 0.3
    trend_str = ', with savings growing with altitude';
elseif r(1, 2) < -0.3
    trend_str = ', with savings growing at lower altitudes';
else
    trend_str = '';
end

abs_sentence = sprintf( ...
    ['Walker-Delta synthesis reduces the minimum satellite count by %.0f\\%%--%.0f\\%% ' ...
     'relative to the Walker-Star \\ac{SOC} baseline across the %d--%d~km altitude range%s.'], ...
    pct_min, pct_max, alt_min, alt_max, trend_str);

fprintf('\n===== Abstract / Conclusion Sentence =====\n');
fprintf('%s\n', abs_sentence);

%% ---- Save .tex file ----------------------------------------------------
tex_file = fullfile(sweep_folder, 'walker_comparison_table.tex');
fid = fopen(tex_file, 'w');
if fid == -1
    warning('Could not write .tex file to %s', tex_file);
else
    for k = 1:numel(tex_lines)
        fprintf(fid, '%s\n', tex_lines{k});
    end
    fprintf(fid, '\n%%%% Abstract / conclusion sentence:\n');
    fprintf(fid, '%% %s\n', abs_sentence);
    fclose(fid);
    fprintf('\nLaTeX table saved to:\n  %s\n', tex_file);
end

%% ---- Multi-stage filtering table ---------------------------------------
print_multistage_table(loaded, sweep_folder, workspace_root);

end

%% ========================================================================
%  Local helpers
%% ========================================================================

function print_multistage_table(loaded, sweep_folder, workspace_root)
% Aggregate multi-stage filtering stats from per-altitude gridsearch runs.

if ~isfield(loaded, 'gridsearch_dirs')
    fprintf('\nNote: gridsearch_dirs not in sweep mat. Re-run sweep to generate multistage table.\n');
    fprintf('(Field is written by the updated Find_valid_constellations.m.)\n');
    return;
end

gridsearch_dirs = loaded.gridsearch_dirs;
total_n1 = 0; total_n2 = 0; total_n3 = 0;
total_wall = 0; total_cands = 0;
total_n_full = 0; total_n_filtered = 0;
sum_t1 = 0; sum_t2 = 0; sum_t3 = 0;
n_workers_rep = 0; n_runs_loaded = 0;

for k = 1:numel(gridsearch_dirs)
    gs_dir = gridsearch_dirs{k};
    if isempty(gs_dir), continue; end
    if gs_dir(1) ~= '/' && ~(numel(gs_dir) > 1 && gs_dir(2) == ':')
        gs_dir = fullfile(workspace_root, gs_dir);
    end
    pd_file = fullfile(gs_dir, 'plot_data.mat');
    if ~isfile(pd_file), continue; end
    tmp = load(pd_file, 'plot_data');
    pd  = tmp.plot_data;
    if ~isfield(pd, 'n_stage1'), continue; end

    total_n1         = total_n1         + pd.n_stage1;
    total_n2         = total_n2         + pd.n_stage2;
    total_n3         = total_n3         + pd.n_stage3;
    total_wall       = total_wall       + pd.wall_time_s;
    total_cands      = total_cands      + pd.candidates_found;
    total_n_full     = total_n_full     + pd.n_full_grid;
    total_n_filtered = total_n_filtered + pd.n_total_grid;
    sum_t1 = sum_t1 + pd.n_stage1 * pd.avg_t1_s;
    sum_t2 = sum_t2 + pd.n_stage2 * pd.avg_t2_s;
    sum_t3 = sum_t3 + pd.n_stage3 * pd.avg_t3_s;
    n_workers_rep = max(n_workers_rep, pd.num_workers);
    n_runs_loaded = n_runs_loaded + 1;
end

if n_runs_loaded == 0
    fprintf('\nNo gridsearch plot_data.mat files with profiling stats found.\n');
    fprintf('Run a new sweep to populate (old runs lack the stats fields).\n');
    return;
end

avg_t1_all = sum_t1 / max(total_n1, 1);
avg_t2_all = sum_t2 / max(total_n2, 1);
avg_t3_all = sum_t3 / max(total_n3, 1);
frac_eval  = 100 * total_n1 / max(total_n_filtered, 1);

fprintf('\n===== Multi-Stage Filtering Summary (%d altitude runs) =====\n', n_runs_loaded);
fprintf('Full grid (all combos):      %s constellations\n', fmt_num(total_n_full));
fprintf('Filtered grid (>=min_sats):  %s constellations\n', fmt_num(total_n_filtered));
fprintf('Stage 1 evaluated:           %s  (%.1f%% of filtered)\n', fmt_num(total_n1), frac_eval);
fprintf('Stage 2 evaluated:           %s  (%.1f%% of S1)\n', fmt_num(total_n2), 100*total_n2/max(total_n1,1));
fprintf('Stage 3 evaluated:           %s  (%.1f%% of S2)\n', fmt_num(total_n3), 100*total_n3/max(total_n2,1));
fprintf('Valid candidates:            %d\n', total_cands);
fprintf('Total wall time:             %.0f s\n', total_wall);
fprintf('Workers:                     %d\n', n_workers_rep);
fprintf('Avg time Stage 1:            %.3f s\n', avg_t1_all);
fprintf('Avg time Stage 2:            %.2f s\n', avg_t2_all);
fprintf('Avg time Stage 3:            %.1f s\n', avg_t3_all);

% Compute "% of previous row" for each stage/filter step
pct_filtered_of_full = 100 * total_n_filtered / max(total_n_full, 1);
pct_s1_of_filtered   = 100 * total_n1         / max(total_n_filtered, 1);
pct_s2_of_s1         = 100 * total_n2         / max(total_n1, 1);
pct_s3_of_s2         = 100 * total_n3         / max(total_n2, 1);
wall_hours           = total_wall / 3600;

ms_tex = {};
ms_tex{end+1} = sprintf('%% Multi-stage filtering table (%d altitude runs, aggregated)', n_runs_loaded);
ms_tex{end+1} = '\begin{table}[htbp]';
ms_tex{end+1} = '    \centering';
ms_tex{end+1} = '    \caption{Multi-Stage Filtering Performance}';
ms_tex{end+1} = '    \label{tab:multistage_results}';
ms_tex{end+1} = '    \begin{tabular}{lccc}';
ms_tex{end+1} = '        \hline';
ms_tex{end+1} = '        \textbf{Stage} & \textbf{Constellations} & \textbf{\% of Prev.} & \textbf{Avg.\ Time (s)} \\ \hline';
ms_tex{end+1} = sprintf('        Full Grid    & %s & -- & -- \\\\',         fmt_num(total_n_full));
ms_tex{end+1} = sprintf('        Filtered Grid & %s & %.1f%% & -- \\\\',  fmt_num(total_n_filtered), pct_filtered_of_full);
ms_tex{end+1} = sprintf('        1 & %s & %.1f\\%% & %.3f \\\\',           fmt_num(total_n1), pct_s1_of_filtered, avg_t1_all);
ms_tex{end+1} = sprintf('        2 & %s & %.1f\\%% & %.2f \\\\',           fmt_num(total_n2), pct_s2_of_s1, avg_t2_all);
ms_tex{end+1} = sprintf('        3 & %s & %.1f\\%% & %.1f \\\\ \\hline',   fmt_num(total_n3), pct_s3_of_s2, avg_t3_all);
ms_tex{end+1} = sprintf('        \\multicolumn{4}{l}{\\small \\textbf{Hardware:} i9-13900K, %d parallel workers.} \\\\', n_workers_rep);
ms_tex{end+1} = sprintf( ...
    ['        \\multicolumn{4}{l}{\\small \\textbf{Total Wall Time:} %.2f~h' ...
     ' \\hspace{0.5cm} \\textbf{Valid Constellations:} %d} \\\\'], wall_hours, total_cands);
ms_tex{end+1} = '        \hline';
ms_tex{end+1} = '    \end{tabular}';
ms_tex{end+1} = '\end{table}';

ms_tex_file = fullfile(sweep_folder, 'multistage_table.tex');
fid2 = fopen(ms_tex_file, 'w');
if fid2 ~= -1
    for k2 = 1:numel(ms_tex)
        fprintf(fid2, '%s\n', ms_tex{k2});
    end
    fclose(fid2);
    fprintf('\nMultistage LaTeX table saved to:\n  %s\n', ms_tex_file);
end
end


function generate_parameter_sweep_plots(heights_km, best_delta_sats, sweep_folder)
% Scatter plots showing how each optimal Walker-Delta parameter distributes
% across orbital altitudes.  Each point = one altitude's best constellation.
% Colorbar encodes orbital altitude throughout.

delta_planes  = [best_delta_sats.Num_planes]';
delta_satspp  = [best_delta_sats.Sats_per_plane]';
delta_phasing = [best_delta_sats.Phasing]';
delta_inc     = [best_delta_sats.Inclination]';

% Convert raw phasing integer F to degrees (consistent with plot_gridsearch)
phasing_deg = (delta_phasing ./ delta_planes) * 360;

alt_clim = [min(heights_km), max(heights_km)];
sz       = 60;
jit      = @(v, s) v + (rand(size(v)) - 0.5) * 2 * s;   % jitter for discrete axes

set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

%% ---- Phasing -----------------------------------------------------------
f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 400]);
scatter(phasing_deg, delta_planes, sz, heights_km, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
colormap(parula); clim(alt_clim);
xlim([0 360]);
xlabel('Phasing (°)', 'FontWeight', 'bold');
ylabel('Number of Planes', 'FontWeight', 'bold');
cb = colorbar; cb.Label.String = 'Orbital Altitude (km)';
title('Optimal Walker-Delta Phasing Across Altitudes');
set(gca, 'FontSize', 14); grid on;
exportgraphics(f1, fullfile(sweep_folder, 'Delta_Phasing_vs_Altitude.png'), 'Resolution', 300);
close(f1);

%% ---- Inclination -------------------------------------------------------
f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 400]);
scatter(delta_inc, heights_km, sz, heights_km, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
colormap(parula); clim(alt_clim);
xlabel('Inclination (°)', 'FontWeight', 'bold');
ylabel('Orbital Altitude (km)', 'FontWeight', 'bold');
title('Optimal Walker-Delta Inclination Across Altitudes');
set(gca, 'FontSize', 14); grid on;
exportgraphics(f2, fullfile(sweep_folder, 'Delta_Inclination_vs_Altitude.png'), 'Resolution', 300);
close(f2);

%% ---- Number of planes --------------------------------------------------
f3 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 400]);
scatter(delta_planes, heights_km, sz, heights_km, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
colormap(parula); clim(alt_clim);
xlabel('Number of Planes', 'FontWeight', 'bold');
ylabel('Orbital Altitude (km)', 'FontWeight', 'bold');
title('Optimal Walker-Delta Planes Across Altitudes');
set(gca, 'FontSize', 14); grid on;
exportgraphics(f3, fullfile(sweep_folder, 'Delta_Planes_vs_Altitude.png'), 'Resolution', 300);
close(f3);

%% ---- Architecture map (sats/plane vs planes) ---------------------------
f4 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 400]);
scatter(delta_planes, delta_satspp, sz, heights_km, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
colormap(parula); clim(alt_clim);
cb = colorbar; cb.Label.String = 'Orbital Altitude (km)';
xlabel('Number of Planes', 'FontWeight', 'bold');
ylabel('Satellites per Plane', 'FontWeight', 'bold');
title('Optimal Walker-Delta Architecture Across Altitudes');
set(gca, 'FontSize', 14); grid on;
exportgraphics(f4, fullfile(sweep_folder, 'Delta_Architecture_Map.png'), 'Resolution', 300);
close(f4);

fprintf('Parameter distribution plots saved to: %s\n', sweep_folder);
end


function s = fmt_num(n)
% Format integer n with thousands separators (e.g. 17791 -> '17,791').
s = num2str(n);
k = length(s) - 3;
while k > 0
    s = [s(1:k) ',' s(k+1:end)];
    k = k - 3;
end
end
