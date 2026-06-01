function plot_constellation_comparison()
% PLOT_CONSTELLATION_COMPARISON
%   Side-by-side 3D globe visualisation of the two optimal 1000 km
%   constellations synthesised in the thesis:
%
%     Walker Delta  – 4 planes × 14 sats, i = 76°  (regional high-inclination)
%     Walker Star   – 5 planes × 14 sats, i = 87°  (near-polar, cross-seam)
%
%   Both constellations serve the same target coverage area (Arctic/Nordic
%   region) but achieve it in fundamentally different ways, making them an
%   ideal pair to illustrate in the report.
%
%   Saves two PNG images to plotting_scripts/figures/constellations/:
%     constellation_walker_delta_1000km.png
%     constellation_walker_star_1000km.png
% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "plot_constellation_comparison"

%% ── Path setup ──────────────────────────────────────────────────────────
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'functions'));  % bootstrap so path_setup is found
path_setup();   % adds project root + all functions/ subfolders to MATLAB path

script_dir  = fileparts(mfilename('fullpath'));
out_dir     = fullfile(script_dir, 'figures', 'constellations');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ── Load optimal constellations ─────────────────────────────────────────
project_root = fileparts(fileparts(script_dir));
opt          = load(fullfile(project_root, 'optimal_constellations.mat'));
target_h_km  = 700;
idx          = find(opt.heights_km == target_h_km, 1, 'first');
if isempty(idx)
    error('No entry for %d km in optimal_constellations.mat', target_h_km);
end

delta_raw = opt.best_delta_sats(idx, :);
star_raw  = opt.star_sats(idx);

fprintf('Walker Delta  – planes: %d  sats/plane: %d  inc: %g°  total: %d\n', ...
    delta_raw.Num_planes, delta_raw.Sats_per_plane, delta_raw.Inclination, delta_raw.Total_sats);
fprintf('Walker Star   – planes: %d  sats/plane: %d  inc: %g°  total: %d\n', ...
    star_raw.Num_planes, star_raw.Sats_per_plane, star_raw.Inclination, star_raw.Total_sats);

%% ── Build full Cfg structs ───────────────────────────────────────────────
% get_cfg resolves optimal_constellations.mat relative to pwd, so we must
% temporarily be in the project root when calling it.
prev_dir = pwd();
cd(project_root);

Cfg_delta              = get_cfg(target_h_km, 'walkerdelta', 'medium', 'short');
Cfg_delta.Num_planes   = delta_raw.Num_planes;
Cfg_delta.Sats_per_plane = delta_raw.Sats_per_plane;
Cfg_delta.Total_sats   = delta_raw.Num_planes * delta_raw.Sats_per_plane;
Cfg_delta.Inclination  = delta_raw.Inclination;
Cfg_delta.Phasing      = delta_raw.Phasing;
Cfg_delta.WalkerStar   = false;

Cfg_star               = get_cfg(target_h_km, 'walkerstar', 'medium', 'short');
Cfg_star.Num_planes    = star_raw.Num_planes;
Cfg_star.Sats_per_plane = star_raw.Sats_per_plane;
Cfg_star.Total_sats    = star_raw.Num_planes * star_raw.Sats_per_plane;
Cfg_star.Inclination   = star_raw.Inclination;
Cfg_star.Phasing       = star_raw.Phasing;
Cfg_star.WalkerStar    = true;

cd(prev_dir);   % restore working directory

%% ── Render and save ──────────────────────────────────────────────────────
constellations = {Cfg_delta, Cfg_star};
labels         = {'walker_delta', 'walker_star'};
titles         = { ...
    sprintf('Walker Delta  |  %d planes × %d sats  |  i = %g°  |  %d km', ...
        Cfg_delta.Num_planes, Cfg_delta.Sats_per_plane, Cfg_delta.Inclination, target_h_km), ...
    sprintf('Walker Star   |  %d planes × %d sats  |  i = %g°  |  %d km', ...
        Cfg_star.Num_planes,  Cfg_star.Sats_per_plane,  Cfg_star.Inclination,  target_h_km)};

for k = 1:2
    Cfg      = constellations{k};
    % Pass empty UE arrays so no ground stations are added to the scene
    Cfg.Flat_UE_array.Lats = [];
    Cfg.Flat_UE_array.Lons = [];
    filename = fullfile(out_dir, sprintf('constellation_%s_%dkm.png', labels{k}, target_h_km));

    fprintf('\n[%d/2] Rendering %s ...\n', k, titles{k});

    show_constellation(Cfg, ...
        false, ...          % show_interactive – headless render
        true,  ...          % save_fig
        out_dir, ...        % out_dir
        false);             % show_details

    % show_constellation saves with its own naming convention; rename to ours
    generated_name = fullfile(out_dir, ...
        sprintf('Constellation_3D_%d_Planes.png', Cfg.Num_planes));
    if isfile(generated_name)
        movefile(generated_name, filename);
        fprintf('    Saved → %s\n', filename);
    else
        fprintf('    [!] Expected file not found: %s\n', generated_name);
    end
end

fprintf('\nDone. Images saved to:\n  %s\n', out_dir);
end
