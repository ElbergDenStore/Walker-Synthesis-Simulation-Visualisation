function out_dir = run_altitude_sweep(Master_config, heights_km, resume_dir)
% RUN_ALTITUDE_SWEEP  Gridsearch the minimal valid constellation per altitude.
%
%   out_dir = run_altitude_sweep(Master_config, heights_km, resume_dir)
%
%   Shared engine behind numerical_walker_synthesis.m. For each orbital
%   altitude (descending) it:
%     1. computes an analytical Walker-Star reference (star_sats),
%     2. runs gridsearch() to find the minimal valid constellation,
%     3. checkpoints results after every altitude so the sweep is resumable.
%
%   Inputs
%     Master_config - struct from synthesis_preset() (+ any overrides). Must
%                     contain field WalkerStar (logical) and folder_prefix.
%     heights_km    - vector of altitudes [km]. Scalar = single-altitude run.
%     resume_dir    - (optional) path to a previous sweep folder to continue.
%                     '' or omitted starts a fresh run.
%
%   Output
%     out_dir       - the sweep output folder under simulation_output/.
%
%   Saved variables use the legacy names (best_delta_sats / best_star_sats,
%   all_delta_sats / all_star_sats) so get_cfg.m and plot_sweep_*.m keep working.

    if nargin < 3 || isempty(resume_dir); resume_dir = ''; end

    is_star      = logical(Master_config.WalkerStar);
    best_field   = ternary(is_star, 'best_star_sats', 'best_delta_sats');
    all_field    = ternary(is_star, 'all_star_sats',  'all_delta_sats');

    heights_km = sort(heights_km, 'descend');

    start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    fprintf('=======================================================\n');
    fprintf('STARTING %s SWEEP AT: %s\n', ...
        ternary(is_star, 'WALKER STAR', 'WALKER DELTA'), char(start_time));
    fprintf('=======================================================\n');

    % --- Set up (or resume) the output folder + checkpoint state ---
    if ~isempty(resume_dir) && exist(resume_dir, 'dir')
        out_dir = resume_dir;
        checkpoint_file = fullfile(out_dir, 'Master_Altitude_Sweep_Results.mat');
        if ~exist(checkpoint_file, 'file')
            error('resume_dir specified but no checkpoint .mat found in: %s', out_dir);
        end
        ck         = load(checkpoint_file);
        best_sats  = ck.(best_field);
        star_sats  = ck.star_sats;

        % altitude -> gridsearch_dir map survives changes to the sweep range
        n_ck = min(numel(ck.heights_km), numel(ck.gridsearch_dirs));
        ck_dir_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
        for ck_i = 1:n_ck
            if ~isempty(ck.gridsearch_dirs{ck_i})
                ck_dir_map(ck.heights_km(ck_i)) = ck.gridsearch_dirs{ck_i};
            end
        end

        gridsearch_dirs = cell(numel(heights_km), 1);
        all_sats        = cell(numel(heights_km), 1);
        for ck_i = 1:numel(heights_km)
            h = heights_km(ck_i);
            if isKey(ck_dir_map, h)
                gridsearch_dirs{ck_i} = ck_dir_map(h);
            end
        end

        completed_heights = heights_km(~cellfun(@isempty, gridsearch_dirs));
        fprintf('Resuming from: %s\n', out_dir);
        fprintf('Already completed %d altitudes. Skipping: %s km\n', ...
            numel(completed_heights), num2str(completed_heights));
    else
        date_str_start  = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name     = sprintf('%s_%s', char(Master_config.folder_prefix), date_str_start);
        out_dir         = fullfile('simulation_output', folder_name);
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end
        star_sats         = [];
        best_sats         = table();
        all_sats          = cell(numel(heights_km), 1);
        gridsearch_dirs   = cell(numel(heights_km), 1);
        completed_heights = [];
    end
    gridsearch_base = fullfile(out_dir, 'gridsearch_runs');

    % --- Altitude loop (descending so min_sats hint propagates downward) ---
    for i = 1:numel(heights_km)

        % Analytical Walker-Star reference for this altitude
        minimum_lat_deg = min(Master_config.Lat_range_deg);
        [Num_planes, Sats_per_plane, Total_sats] = calculate_walker_star( ...
            heights_km(i), minimum_lat_deg, Master_config.Min_elevation_UE);
        star_sats(i).Orbit_height   = heights_km(i); %#ok<AGROW>
        star_sats(i).Total_sats     = Total_sats;
        star_sats(i).Num_planes     = Num_planes;
        star_sats(i).Phasing        = Num_planes / 2;
        star_sats(i).Inclination    = 90;
        star_sats(i).Sats_per_plane = Sats_per_plane;

        if ismember(heights_km(i), completed_heights)
            fprintf('Skipping already-completed altitude %d km\n', heights_km(i));
            continue;
        end

        % Constrain the search using the best result from the previous (higher) altitude
        if height(best_sats) > 0
            min_sats = best_sats.Total_sats(end);
        else
            min_sats = 0;
        end

        [best_params, all_sats{i}, gridsearch_dirs{i}] = gridsearch( ...
            Master_config, heights_km(i), min_sats, gridsearch_base);

        best_sats = [best_sats; best_params(1, :)]; %#ok<AGROW>

        % --- CHECKPOINT after every completed altitude ---
        save_checkpoint(out_dir, heights_km, star_sats, best_sats, all_sats, ...
            Master_config, gridsearch_dirs, best_field, all_field);
        fprintf('[Checkpoint] Saved after altitude %d km\n', heights_km(i));
    end

    % --- Final save ---
    save_checkpoint(out_dir, heights_km, star_sats, best_sats, all_sats, ...
        Master_config, gridsearch_dirs, best_field, all_field);
    fprintf('Results saved to: %s\n', out_dir);

    % Walker Delta sweeps feed get_cfg.m via optimal_constellations.mat
    if ~is_star
        best_delta_sats = best_sats; %#ok<NASGU>
        save('optimal_constellations.mat', 'heights_km', 'star_sats', 'best_delta_sats');
        fprintf('Optimal constellations saved to optimal_constellations.mat\n');
    end

    end_time     = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    elapsed_time = end_time - start_time;
    fprintf('\n=======================================================\n');
    fprintf('SWEEP FINISHED AT: %s\n', char(end_time));
    fprintf('TOTAL ELAPSED TIME: %s\n', char(elapsed_time));
    fprintf('=======================================================\n');
end

% -------------------------------------------------------------------------
function save_checkpoint(out_dir, heights_km, star_sats, best_sats, all_sats, ...
        Master_config, gridsearch_dirs, best_field, all_field) %#ok<INUSD>
    S = struct();
    S.heights_km        = heights_km;
    S.star_sats         = star_sats;
    S.Master_config     = Master_config;
    S.gridsearch_dirs   = gridsearch_dirs;
    S.(best_field)      = best_sats;
    S.(all_field)       = all_sats;
    save(fullfile(out_dir, 'Master_Altitude_Sweep_Results.mat'), '-struct', 'S');
end

% -------------------------------------------------------------------------
function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end
