%% compare_constellation_population_coverage.m
% Compares Walker-Star, global Walker-Delta, and SPLIT Walker-Delta
% (low-lat 0-55 + high-lat 55-85, summed) constellations on a per-latitude
% basis, weighted by global population density (WorldPop 2020, 1 km).
%
% Self-contained: reads pre-computed constellation tables from data/*.mat.
% Adds repo root + functions/ to the path so constellation_simulator
% and helpers are available regardless of caller cwd.
%
% Outputs:
%   figures/Coverage_vs_Population.png
%   tables/pop_weighted_coverage_table.tex
%   results_avg_visible.mat   (raw per-scenario per-altitude visibility curves)
%
% Author: Jacob Elberg Nielsen, May 2026

clearvars; close all; clc;

%% ===== Paths =====
this_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(this_dir, 'data');
fig_dir  = fullfile(this_dir, 'figures');
tab_dir  = fullfile(this_dir, 'tables');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
if ~exist(tab_dir,'dir'), mkdir(tab_dir); end

% Repo root is two folders up (plotting_scripts/coverage_population_comparison)
repo_root = fileparts(fileparts(this_dir));
addpath(repo_root);
addpath(fullfile(repo_root, 'functions'));
path_setup();   % adds data/ etc. to MATLAB path so readgeoraster finds the tif

%% ===== Run config =====
altitudes_km   = [500, 700, 1000, 1200];
lat_step_deg   = 0.25;        % probe-grid latitude resolution
probe_lons_deg = [-135, -45, 45, 135];   % 4 longitudes; results averaged
sim_hours      = 12;
sample_time_s  = 60;
min_elev_deg   = 20;
% Serial execution avoids broadcast-variable memory blowup.
% constellation_simulator with use_parallel=true sends sat_pos_ecef
% to every worker; after 3-4 large sims the workers run out of heap.
% Each per-UE iteration is fast vectorised geometry, so serial is fine.
use_parallel   = false;

%% ===== Figure style =====
FIG.dpi        = 300;
FIG.font_size  = 14;
FIG.title_size = 14;
FIG.lw         = 1.8;
FIG.fig_size   = [100 100 600 400];   % each individual figure
set(groot,'DefaultAxesFontSize',FIG.font_size);
set(groot,'DefaultTextFontSize',FIG.font_size);
set(groot,'DefaultLineLineWidth',FIG.lw);

%% ===== Load global population GeoTIFF =====
fprintf('Loading global population GeoTIFF...\n');
tif_path = which('ppp_2020_1km_Aggregated.tif');
if isempty(tif_path)
    error('ppp_2020_1km_Aggregated.tif not found on path. Ensure functions/data is on the path.');
end
[pop_data, R] = readgeoraster(tif_path);
pop_data = double(pop_data);
pop_data(pop_data < 0 | isnan(pop_data)) = 0;
tif_lat_lim = R.LatitudeLimits;
fprintf('  TIF latitude coverage: [%.2f, %.2f] deg\n', tif_lat_lim(1), tif_lat_lim(2));

% Probe grid: always symmetric ±90° so north/south mirrors always exist.
% Population bins outside the TIF latitude range just sum to zero (handled below).
probe_lats = (-90:lat_step_deg:90).';
nProbeLats = numel(probe_lats);
nProbeLons = numel(probe_lons_deg);

% Compute population per pixel row (sum over all longitudes = full globe)
nRows = size(pop_data, 1);
center_col_vec = ones(nRows,1) * round(size(pop_data,2)/2);
[row_lats, ~] = intrinsicToGeographic(R, center_col_vec, (1:nRows).');
pop_per_row = sum(pop_data, 2);

% Bin into probe_lats grid
pop_per_probe_lat = zeros(size(probe_lats));
for i = 1:nProbeLats
    edges_lo = probe_lats(i) - lat_step_deg/2;
    edges_hi = probe_lats(i) + lat_step_deg/2;
    mask = row_lats >= edges_lo & row_lats < edges_hi;
    pop_per_probe_lat(i) = sum(pop_per_row(mask));
end
fprintf('  Total binned population (full globe): %.2e\n', sum(pop_per_probe_lat));

%% ===== Load constellation parameter tables (self-contained data/) =====
fprintf('Loading constellation parameter files...\n');
src_star    = load(fullfile(data_dir,'star_global.mat'));
src_dglobal = load(fullfile(data_dir,'delta_global.mat'));
src_dlo     = load(fullfile(data_dir,'delta_lowlat.mat'));
src_dhi     = load(fullfile(data_dir,'delta_highlat.mat'));

scenarios = struct( ...
    'name',   {'Walker-Star (global)',   'Walker-Delta (global)', 'Walker-Delta (0-55)', 'Walker-Delta (55-90)'}, ...
    'key',    {'star',                   'delta_global',          'delta_lowlat',        'delta_highlat'}, ...
    'is_star',{true,                     false,                   false,                 false}, ...
    'src',    {src_star,                 src_dglobal,             src_dlo,               src_dhi});
nSc = numel(scenarios);

%% ===== Probe UE grid (flat lat,lon list) =====
[LAT, LON]    = ndgrid(probe_lats, probe_lons_deg);
probe_UE_lats = LAT(:);
probe_UE_lons = LON(:);

%% ===== Run simulations (checkpoint/resume) =====
nAlt            = numel(altitudes_km);
avg_visible_mat = nan(nProbeLats, nSc, nAlt);
total_sats_mat  = nan(nSc, nAlt);
checkpoint_file = fullfile(this_dir, 'checkpoint_avg_visible.mat');

% --- Resume from checkpoint if present ---
if isfile(checkpoint_file)
    ck = load(checkpoint_file, 'avg_visible_mat', 'total_sats_mat');
    if size(ck.avg_visible_mat, 1) == nProbeLats
        % Grid matches — direct load.
        avg_visible_mat = ck.avg_visible_mat;
        total_sats_mat  = ck.total_sats_mat;
        fprintf('[RESUME] Loaded checkpoint. Skipping already-completed scenarios.\n');
    else
        % Grid size mismatch (old checkpoint used TIF-restricted probe_lats).
        % Migrate: copy existing latitudes, mirror north→south for the gap.
        old_size = size(ck.avg_visible_mat, 1);
        old_lats = (90 - (old_size-1)*lat_step_deg : lat_step_deg : 90).';
        old_mat  = ck.avg_visible_mat;
        fprintf('[MIGRATE] Checkpoint has %d lats (%.1f° to %.1f°); extending to %d (-90° to 90°) by symmetry.\n', ...
            old_size, old_lats(1), old_lats(end), nProbeLats);
        new_avg = nan(nProbeLats, nSc, nAlt);
        for ii = 1:nProbeLats
            lat = probe_lats(ii);
            d1  = abs(old_lats - lat);
            if min(d1) < lat_step_deg/2
                new_avg(ii,:,:) = old_mat(d1 == min(d1),:,:);
            else
                d2 = abs(old_lats - (-lat));
                if min(d2) < lat_step_deg/2
                    new_avg(ii,:,:) = old_mat(d2 == min(d2),:,:);
                end
            end
        end
        avg_visible_mat = new_avg;
        total_sats_mat  = ck.total_sats_mat;
        fprintf('[MIGRATE] Done.\n');
    end
    clear ck
end

for ai = 1:nAlt
    h_km = altitudes_km(ai);
    fprintf('\n==== Altitude %d km ====\n', h_km);

    for si = 1:nSc
        % Skip if already done (checkpoint)
        if ~isnan(avg_visible_mat(1, si, ai))
            fprintf('  [skip] %s @ %d km (already done)\n', scenarios(si).name, h_km);
            continue
        end

        sc_def = scenarios(si);
        fprintf('  [%d/%d] Scenario: %s  ...\n', (ai-1)*nSc+si, nAlt*nSc, sc_def.name);
        t_sim = tic;

        Cfg = build_cfg_for_scenario(sc_def, h_km, sample_time_s, min_elev_deg, sim_hours);
        Cfg.Flat_UE_array.Lats = probe_UE_lats;
        Cfg.Flat_UE_array.Lons = probe_UE_lons;
        Cfg.NumUEs = numel(probe_UE_lats);

        m = constellation_simulator(Cfg, use_parallel, false);

        % Extract ONLY Num_visible, then immediately free the full metrics struct.
        % metrics.UEs contains all time series (Time, Range, Elevation, Azimuth,
        % SatID) for every UE — hundreds of MB — which accumulates across sims.
        num_vis    = m.Num_visible;  % [NumUEs x nT] — this is all we need
        n_sats_now = Cfg.Total_sats; % capture before clear
        clear m Cfg

        avg_per_ue = mean(num_vis, 2);
        avg_per_ue = reshape(avg_per_ue, nProbeLats, nProbeLons);
        avg_visible_mat(:, si, ai) = mean(avg_per_ue, 2);
        total_sats_mat(si, ai)     = n_sats_now;
        clear num_vis avg_per_ue n_sats_now

        fprintf('    done in %.1f s\n', toc(t_sim));

        % Checkpoint after every scenario
        save(checkpoint_file, 'avg_visible_mat', 'total_sats_mat', 'probe_lats');
    end
end
fprintf('\nAll simulations complete. Checkpoint saved.\n');

%% ===== Save raw results =====
save(fullfile(this_dir,'results_avg_visible.mat'), ...
    'probe_lats','probe_lons_deg','altitudes_km','scenarios', ...
    'avg_visible_mat','total_sats_mat','pop_per_probe_lat','tif_lat_lim', ...
    'sim_hours','sample_time_s','min_elev_deg');

%% ===== Derived: split delta = low-lat + high-lat =====
idx_st = find(strcmp({scenarios.key},'star'));
idx_gl = find(strcmp({scenarios.key},'delta_global'));
idx_lo = find(strcmp({scenarios.key},'delta_lowlat'));
idx_hi = find(strcmp({scenarios.key},'delta_highlat'));

split_visible    = squeeze(avg_visible_mat(:,idx_lo,:) + avg_visible_mat(:,idx_hi,:));   % [nLats x nAlt]
split_total_sats = total_sats_mat(idx_lo,:) + total_sats_mat(idx_hi,:);                  % [1 x nAlt]

%% ===== FIGURE: 4 separate figures, sine-latitude equal-area x-axis =====
col_star  = [0.00 0.45 0.74];
col_glob  = [0.85 0.33 0.10];
col_split = [0.47 0.67 0.19];

pop_smooth = smoothdata(pop_per_probe_lat, 'gaussian', 11);
x_sine     = sind(probe_lats);

% Symmetric ticks: [-90,-60,-30,0,30,60,90] map to [-1,-0.866,-0.5,0,0.5,0.866,1]
% which are evenly distributed on the sine axis.
tick_lats = [-90, -60, -30, 0, 30, 60, 90];
tick_pos  = sind(tick_lats);
tick_lbls = arrayfun(@(d) [num2str(d) char(176)], tick_lats, 'UniformOutput', false);

% Symmetrize visibility curves: constellations are equatorially symmetric;
% any difference between vis(+phi) and vis(-phi) is simulation noise.
% For each probe latitude, average with its mirror latitude where both exist.
lat_step_val = probe_lats(2) - probe_lats(1);
mirror_idx   = round((-probe_lats - probe_lats(1)) / lat_step_val) + 1;
valid_mirror = mirror_idx >= 1 & mirror_idx <= nProbeLats;

for ai = 1:nAlt
    vis_star  = symmetrize_vis(avg_visible_mat(:,idx_st,ai), mirror_idx, valid_mirror);
    vis_glob  = symmetrize_vis(avg_visible_mat(:,idx_gl,ai), mirror_idx, valid_mirror);
    vis_split = symmetrize_vis(split_visible(:,ai),          mirror_idx, valid_mirror);

    vis_max = max([vis_star; vis_glob; vis_split]);
    if isnan(vis_max) || vis_max <= 0; vis_max = 1; end
    y_top = vis_max * 1.18;

    fig_i = figure('Color','w','Position', FIG.fig_size);
    hold on; grid on; box on;

    % Population density background (area-preserving on sine axis: divide by cos(lat))
    cos_lats       = max(cosd(probe_lats), 0.01);
    pop_per_sinlat = pop_smooth ./ cos_lats;
    pop_norm       = pop_per_sinlat / max(pop_per_sinlat);
    area(x_sine, pop_norm * y_top * 0.70, ...
        'FaceColor',[0.78 0.78 0.78],'FaceAlpha',0.55,'EdgeColor','none', ...
        'DisplayName', 'Pop. density (normalised, area-preserving)');

    % Visibility curves (symmetrized)
    plot(x_sine, vis_star,  '-',  'Color', col_star,  'LineWidth', FIG.lw, ...
        'DisplayName', ['Walker-Star  (T=' num2str(total_sats_mat(idx_st,ai)) ')']);
    plot(x_sine, vis_glob,  '-', 'Color', col_glob,  'LineWidth', FIG.lw, ...
        'DisplayName', ['Walker-\Delta global  (T=' num2str(total_sats_mat(idx_gl,ai)) ')']);
    plot(x_sine, vis_split, '-', 'Color', col_split, 'LineWidth', FIG.lw, ...
        'DisplayName', ['Walker-\Delta split  (T=' num2str(split_total_sats(ai)) ')']);

    ylabel('Avg. satellites in view');
    xlabel('Latitude');
    xlim([-1 1]);
    ylim([0, y_top]);
    xticks(tick_pos);
    xticklabels(tick_lbls);
    title(sprintf('h = %d km', altitudes_km(ai)));
    legend('Location','northwest','FontSize', FIG.font_size - 3);

    fig_out = fullfile(fig_dir, sprintf('Coverage_vs_Population_%dkm.png', altitudes_km(ai)));
    exportgraphics(fig_i, fig_out, 'Resolution', FIG.dpi);
    fprintf('  Figure saved: %s\n', fig_out);
end
fprintf('\nAll figures saved.\n');

%% ===== TABLE: constellation parameters (T, P, i) =====
% Extract parameters directly from the source .mat files at each altitude.
star_Tv = nan(nAlt,1); star_Pv = nan(nAlt,1); star_Iv = nan(nAlt,1);
gl_Tv   = nan(nAlt,1); gl_Pv   = nan(nAlt,1); gl_Iv   = nan(nAlt,1);
lo_Tv   = nan(nAlt,1); lo_Pv   = nan(nAlt,1); lo_Iv   = nan(nAlt,1);
hi_Tv   = nan(nAlt,1); hi_Pv   = nan(nAlt,1); hi_Iv   = nan(nAlt,1);

for ai = 1:nAlt
    h = altitudes_km(ai);

    [~, ii] = min(abs(src_star.heights_km(:) - h));
    s = src_star.star_sats(ii);
    star_Tv(ai) = s.Total_sats;  star_Pv(ai) = s.Num_planes;  star_Iv(ai) = s.Inclination;

    [~, ii] = min(abs(src_dglobal.heights_km(:) - h));
    [gl_Tv(ai), gl_Pv(ai), gl_Iv(ai)] = get_delta_params(src_dglobal.best_delta_sats, ii);

    [~, ii] = min(abs(src_dlo.heights_km(:) - h));
    [lo_Tv(ai), lo_Pv(ai), lo_Iv(ai)] = get_delta_params(src_dlo.best_delta_sats, ii);

    [~, ii] = min(abs(src_dhi.heights_km(:) - h));
    [hi_Tv(ai), hi_Pv(ai), hi_Iv(ai)] = get_delta_params(src_dhi.best_delta_sats, ii);
end
split_Tv = lo_Tv + hi_Tv;

row_strs = cell(nAlt, 1);

% Compute max T/P string length per column so inclinations align across rows.
tp_len   = @(T, P) length(sprintf('%d/%d', T, P));
star_max = max(arrayfun(tp_len, star_Tv, star_Pv));
gl_max   = max(arrayfun(tp_len, gl_Tv,   gl_Pv));
lo_max   = max(arrayfun(tp_len, lo_Tv,   lo_Pv));
hi_max   = max(arrayfun(tp_len, hi_Tv,   hi_Pv));

for ai = 1:nAlt
    % Each cell: "T/P\,\,... i°" — thin-space pad so inclinations line up.
    fmt = @(T, P, I, mx) sprintf('%d/%d%s %.0f\\textdegree{}', ...
        T, P, repmat('\,', 1, max(0, mx - tp_len(T,P)) * 2), I);
    row_strs{ai} = sprintf( ...
        '    %d & %s & %s & %d & %s & %s \\\\', ...
        altitudes_km(ai), ...
        fmt(star_Tv(ai), star_Pv(ai), star_Iv(ai), star_max), ...
        fmt(gl_Tv(ai),   gl_Pv(ai),   gl_Iv(ai),   gl_max), ...
        split_Tv(ai), ...
        fmt(lo_Tv(ai),   lo_Pv(ai),   lo_Iv(ai),   lo_max), ...
        fmt(hi_Tv(ai),   hi_Pv(ai),   hi_Iv(ai),   hi_max));
end

tex_lines = [
    {'% Auto-generated by compare_constellation_population_coverage.m'}
    {'\begin{table}[!t]'}
    {'  \caption{Optimal constellation parameters at each orbital altitude. Each entry shows $T$/$P$ (total satellites/planes) and inclination $i$. The composite constellation combines two regional Walker-$\Delta$ constellations (0--55\textdegree{} and 55--90\textdegree{}); $T$ is their combined count.}'}
    {'  \label{tab:constellation_params}'}
    {'  \centering'}
    {'  \renewcommand{\arraystretch}{1.2}'}
    {'  \begin{tabular}{c|l|l|r|l|l}'}
    {'    \hline\hline'}
    {'    & \textbf{Star} & \textbf{Delta} & \multicolumn{3}{c}{\textbf{Composite}} \\'}
    {'    \textbf{Alt.~(km)} & $T$/$P$ $i$ & $T$/$P$ $i$ & $T$ & $0^\circ$--$55^\circ$ & $55^\circ$--$90^\circ$ \\'}
    {'    \hline'}
    row_strs
    {'    \hline\hline'}
    {'  \end{tabular}'}
    {'\end{table}'}
];
tex_path = fullfile(tab_dir, 'constellation_params_table.tex');
fid = fopen(tex_path, 'w');
fprintf(fid, '%s\n', tex_lines{:});
fclose(fid);
fprintf('Table saved: %s\n', tex_path);

fprintf('\nAll done.\n');


%% =================================================================
function Cfg = build_cfg_for_scenario(sc_def, h_km, sample_time_s, min_elev_deg, sim_hours)
% Build a complete Cfg by starting from get_cfg's defaults (which fill all
% link / FRF / etc. fields) and overriding only the constellation params.
    if sc_def.is_star
        Cfg = get_cfg(h_km, "walkerstar");
    else
        Cfg = get_cfg(h_km, "walkerdelta");
    end

    % --- Override scenario constellation from the scenario's source .mat ---
    src = sc_def.src;
    heights = src.heights_km(:);
    [~, idx] = min(abs(heights - h_km));

    if sc_def.is_star
        s = src.star_sats(idx);
        Cfg.Num_planes     = s.Num_planes;
        Cfg.Sats_per_plane = s.Sats_per_plane;
        Cfg.Inclination    = s.Inclination;
        Cfg.WalkerStar     = true;
        Cfg.Total_sats     = Cfg.Num_planes * Cfg.Sats_per_plane;
        Cfg.Phasing        = Cfg.Num_planes / 2;
    else
        t = src.best_delta_sats;
        if istable(t)
            row = t(idx, :);
            Cfg.Num_planes     = row.Num_planes;
            Cfg.Sats_per_plane = row.Sats_per_plane;
            Cfg.Inclination    = row.Inclination;
            Cfg.Phasing        = row.Phasing;
        else
            s = t(idx);
            Cfg.Num_planes     = s.Num_planes;
            Cfg.Sats_per_plane = s.Sats_per_plane;
            Cfg.Inclination    = s.Inclination;
            Cfg.Phasing        = s.Phasing;
        end
        Cfg.WalkerStar = false;
        Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
    end

    % --- Override timing / sim region for this study ---
    Cfg.Orbit_height     = h_km * 1e3;
    Cfg.SampleTime       = sample_time_s;
    Cfg.Min_elevation_UE = min_elev_deg;
    Cfg.StartTime        = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
    Cfg.StopTime         = Cfg.StartTime + hours(sim_hours);
    Cfg.Lat_range_deg    = [-89, 89];   % full-globe probes; affects seam-ratio geometry only
end

%% =================================================================
function vs = symmetrize_vis(v, mirror_idx, valid_mirror)
% Average each latitude bin with its equatorial mirror to remove simulation noise.
% Constellations are symmetric about the equator; any asymmetry is stochastic.
    vs = v;
    for ii = 1:numel(v)
        if valid_mirror(ii)
            vs(ii) = mean([v(ii), v(mirror_idx(ii))]);
        end
    end
end

%% =================================================================
function [T, P, I] = get_delta_params(best_delta_sats, idx)
% Extract T, P, I from a Walker-Delta result that may be a table or struct array.
    if istable(best_delta_sats)
        row = best_delta_sats(idx, :);
        T = row.Total_sats;  P = row.Num_planes;  I = row.Inclination;
    else
        s = best_delta_sats(idx);
        T = s.Total_sats;    P = s.Num_planes;    I = s.Inclination;
    end
end
