function plot_sweep_star(sweep_folder)
% PLOT_SWEEP_STAR  Plot and summarize a Walker Star gridsearch sweep.
%   Compares the analytical Walker-Star formula against the numerical
%   gridsearch result at each altitude.
%
%   Outputs (saved to sweep folder):
%     Star_Analytical_vs_Numerical.png  – scatter comparison plot
%     star_comparison_table.tex          – IEEE LaTeX table
%     multistage_table.tex               – multi-stage filtering stats
%
% Usage:
%   plot_sweep_star()                  – most recent Master_Sweep_WalkerStar_* folder
%   plot_sweep_star('path/to/sweep')   – specific sweep folder

script_dir     = fileparts(mfilename('fullpath'));
workspace_root = fileparts(script_dir);
sim_dir        = fullfile(workspace_root, 'simulation_output');

%% ---- Locate sweep folder -----------------------------------------------
if nargin == 0
    hits = dir(fullfile(sim_dir, 'Master_Sweep_WalkerStar_*'));
    hits = hits([hits.isdir]);
    if isempty(hits)
        error('No Master_Sweep_WalkerStar_* folders found in:\n  %s', sim_dir);
    end
    [~, idx]     = max([hits.datenum]);
    sweep_folder = fullfile(sim_dir, hits(idx).name);
    fprintf('Using most recent Walker Star sweep: %s\n', sweep_folder);
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
if ~isfield(loaded, 'best_star_sats')
    error(['This folder contains a Walker Delta sweep, not a Walker Star sweep.\n' ...
           'Use plot_sweep_delta() instead.']);
end

heights_km     = loaded.heights_km;
star_sats      = loaded.star_sats;       % analytical
best_star_sats = loaded.best_star_sats;  % numerical gridsearch

% Sort ascending for consistent indexing
[heights_km, sort_idx] = sort(heights_km, 'ascend');
star_sats      = star_sats(sort_idx);
best_star_sats = best_star_sats(sort_idx, :);

anal_T = [star_sats.Total_sats]';
num_T  = best_star_sats.Total_sats;

%% ---- Plot ---------------------------------------------------------------
if isfield(loaded, 'Master_config')
    cfg             = loaded.Master_config;
    minimum_lat_deg = min(cfg.Lat_range_deg);
    maximum_lat_deg = max(cfg.Lat_range_deg);
    title_str = sprintf( ...
        'Walker Star: Analytical vs Numerical | Lat: %.1f^{\\circ}--%.1f^{\\circ} | \\epsilon_{min}: %.1f^{\\circ}', ...
        minimum_lat_deg, maximum_lat_deg, cfg.Min_elevation_UE);
else
    title_str = 'Walker Star: Analytical vs Numerical';
end

set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

f1 = figure('Visible', 'off', 'Name', 'Walker Star Comparison', 'Color', 'w', 'Position', [100 100 700 450]);
hold on;
scatter(heights_km, anal_T, 36, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'Analytical');
scatter(heights_km, num_T,  36, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Numerical');
xlabel('Orbit Height (km)', 'FontWeight', 'bold');
ylabel('Total Satellites Required', 'FontWeight', 'bold');
title(title_str);
legend('Location', 'northeast');
grid on; hold off;

plot_path = fullfile(sweep_folder, 'Star_Analytical_vs_Numerical.png');
exportgraphics(f1, plot_path, 'Resolution', 300);
close(f1);
fprintf('Plot saved to: %s\n', plot_path);

%% ---- Console summary ---------------------------------------------------
n          = numel(heights_km);
abs_saving = anal_T - num_T;
pct_saving = 100 * abs_saving ./ anal_T;

fprintf('\n===== Walker-Star Analytical vs Numerical Summary =====\n');
fprintf('%-10s  %-8s  %-8s  %-10s  %-10s\n', ...
    'Alt (km)', 'Anal-T', 'Num-T', 'Abs Save', 'Pct Save (%)');
fprintf('%s\n', repmat('-', 1, 54));
for k = 1:n
    fprintf('%-10d  %-8d  %-8d  %-10d  %-10.1f\n', ...
        heights_km(k), anal_T(k), num_T(k), abs_saving(k), pct_saving(k));
end

fprintf('\nOverall range (all altitudes):\n');
fprintf('  Abs savings : %d – %d satellites\n', min(abs_saving), max(abs_saving));
fprintf('  Pct savings : %.1f%% – %.1f%%\n',    min(pct_saving), max(pct_saving));
if any(pct_saving < 0)
    neg_alts = heights_km(pct_saving < 0);
    fprintf('  NOTE: %d altitude(s) where numerical result exceeds analytical:\n', numel(neg_alts));
    for ki = 1:numel(neg_alts)
        ni = find(heights_km == neg_alts(ki), 1);
        fprintf('    %d km  (Anal=%d, Num=%d) — inclination constraint may differ\n', ...
            neg_alts(ki), anal_T(ni), num_T(ni));
    end
end

%% ---- Representative altitude detail -----------------------------------
rep_alts = [500, 700, 1000, 1200];
rep_idx  = zeros(1, numel(rep_alts));
for k = 1:numel(rep_alts)
    [~, rep_idx(k)] = min(abs(heights_km - rep_alts(k)));
end

fprintf('\n===== Representative Altitude Detail =====\n');
hdr = sprintf('%-6s | %-4s %-4s %-4s %-5s | %-4s %-4s %-4s %-5s %-5s | %-8s %-8s', ...
    'Alt', 'T(A)', 'P(A)', 'S(A)', 'i(A)', ...
    'T(N)', 'P(N)', 'S(N)', 'F(N)', 'i(N)', ...
    'ΔT', 'Save%');
fprintf('%s\n%s\n', hdr, repmat('-', 1, numel(hdr)));
for k = 1:numel(rep_idx)
    ri = rep_idx(k);
    as = star_sats(ri);
    ns = best_star_sats(ri, :);
    fprintf('%-6d | %-4d %-4d %-4d %-5.1f | %-4d %-4d %-4d %-5.2f %-5.2f | %-8d %-8.1f\n', ...
        heights_km(ri), ...
        as.Total_sats, as.Num_planes, as.Sats_per_plane, as.Inclination, ...
        ns.Total_sats, ns.Num_planes, ns.Sats_per_plane, ns.Phasing, ns.Inclination, ...
        as.Total_sats - ns.Total_sats, ...
        100 * (as.Total_sats - ns.Total_sats) / as.Total_sats);
end

%% ---- LaTeX table -------------------------------------------------------
tex_lines = {};
tex_lines{end+1} = '% Auto-generated by plot_sweep_star.m';
tex_lines{end+1} = '% Paste into your IEEE Access manuscript.';
tex_lines{end+1} = '%';
tex_lines{end+1} = '% Column key (A = Analytical, N = Numerical):';
tex_lines{end+1} = '%   T = total satellites, P = planes, S = sats/plane,';
tex_lines{end+1} = '%   F = phasing parameter, i = inclination (degrees)';
tex_lines{end+1} = '';
tex_lines{end+1} = '\begin{table}[!t]';
tex_lines{end+1} = '  \caption{Walker Star: Analytical vs Numerical at Representative Altitudes}';
tex_lines{end+1} = '  \label{tab:star_anal_vs_num}';
tex_lines{end+1} = '  \centering';
tex_lines{end+1} = '  \renewcommand{\arraystretch}{1.2}';
tex_lines{end+1} = '  \begin{tabular}{c|cccc|ccccc|cc}';
tex_lines{end+1} = '    \hline\hline';
tex_lines{end+1} = '    & \multicolumn{4}{c|}{\textbf{Analytical Walker Star}}';
tex_lines{end+1} = '    & \multicolumn{5}{c|}{\textbf{Numerical Walker Star}}';
tex_lines{end+1} = '    & \multicolumn{2}{c}{\textbf{Saving}} \\';
tex_lines{end+1} = '    \textbf{Alt.} & $T$ & $P$ & $S$ & $i$ (°)';
tex_lines{end+1} = '    & $T$ & $P$ & $S$ & $F$ & $i$ (°)';
tex_lines{end+1} = '    & $\Delta T$ & \% \\';
tex_lines{end+1} = '    \textbf{(km)} & & & & & & & & & & & \\';
tex_lines{end+1} = '    \hline';

for k = 1:numel(rep_idx)
    ri  = rep_idx(k);
    as  = star_sats(ri);
    ns  = best_star_sats(ri, :);
    dT  = as.Total_sats - ns.Total_sats;
    pct = 100 * dT / as.Total_sats;
    tex_lines{end+1} = sprintf( ...   %#ok<AGROW>
        '    %d & %d & %d & %d & %.0f & %d & %d & %d & %.2f & %.2f & %d & %.1f\\%% \\\\', ...
        heights_km(ri), ...
        as.Total_sats, as.Num_planes, as.Sats_per_plane, as.Inclination, ...
        ns.Total_sats, ns.Num_planes, ns.Sats_per_plane, ns.Phasing, ns.Inclination, ...
        dT, pct);
end

tex_lines{end+1} = '    \hline\hline';
tex_lines{end+1} = '  \end{tabular}';
tex_lines{end+1} = '\end{table}';

%% ---- Abstract / conclusion sentence ------------------------------------
pct_min  = min(pct_saving);
pct_max  = max(pct_saving);
alt_min  = min(heights_km);
alt_max  = max(heights_km);
[r, ~]   = corrcoef(heights_km, pct_saving);
if r(1, 2) > 0.3
    trend_str = ', with savings growing with altitude';
elseif r(1, 2) < -0.3
    trend_str = ', with savings growing at lower altitudes';
else
    trend_str = '';
end

abs_sentence = sprintf( ...
    ['Numerical Walker-Star optimisation reduces the minimum satellite count by ' ...
     '%.0f\\%%--%.0f\\%% relative to the analytical formula across ' ...
     'the %d--%d~km altitude range%s.'], ...
    pct_min, pct_max, alt_min, alt_max, trend_str);

fprintf('\n===== Abstract / Conclusion Sentence =====\n');
fprintf('%s\n', abs_sentence);

%% ---- Save .tex file ----------------------------------------------------
tex_file = fullfile(sweep_folder, 'star_comparison_table.tex');
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

ms_tex = {};
ms_tex{end+1} = sprintf('%% Multi-stage filtering table (%d altitude runs, aggregated)', n_runs_loaded);
ms_tex{end+1} = '\begin{table}[htbp]';
ms_tex{end+1} = '    \centering';
ms_tex{end+1} = '    \caption{Multi-Stage Constellation Filtering: Configuration and Performance}';
ms_tex{end+1} = '    \label{tab:multistage_results}';
ms_tex{end+1} = '    \begin{tabular}{lcc}';
ms_tex{end+1} = '        \multicolumn{3}{c}{\textbf{Multi Stage Filtering Results}} \\';
ms_tex{end+1} = '        \hline';
ms_tex{end+1} = '        \textbf{Stage} & \textbf{Constellations Evaluated} & \textbf{Avg.\ Time (s)} \\ \hline';
ms_tex{end+1} = sprintf('        1 & %s & %.3f \\\\', fmt_num(total_n1), avg_t1_all);
ms_tex{end+1} = sprintf('        2 & %s & %.2f \\\\', fmt_num(total_n2), avg_t2_all);
ms_tex{end+1} = sprintf('        3 & %s & %.1f \\\\ \\hline', fmt_num(total_n3), avg_t3_all);
ms_tex{end+1} = sprintf('        \\multicolumn{3}{l}{\\small \\textbf{Hardware:} %d parallel workers.} \\\\', n_workers_rep);
ms_tex{end+1} = sprintf( ...
    ['        \\multicolumn{3}{l}{\\small \\textbf{Total Wall Time:} %.0f~s' ...
     ' \\hspace{0.5cm} \\textbf{Valid Candidates:} %d} \\\\'], total_wall, total_cands);
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


function s = fmt_num(n)
% Format integer n with thousands separators (e.g. 17791 -> '17,791').
s = num2str(n);
k = length(s) - 3;
while k > 0
    s = [s(1:k) ',' s(k+1:end)];
    k = k - 3;
end
end
