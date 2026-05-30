% How to run through the night:
% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "Find_valid_constellations" > log.txt
% xvfb is a virtual display to avoid constellation pictures do not crash server
% one ">" overwrites the file
% Read log during run using tail -f log.txt

% Read log afterwards using less
% How to use it: Type less log.txt
% Pro-tips inside less:% 
% Press Space to page down, b to page up.
% Press G to jump immediately to the very bottom (the newest logs).
% Press g to jump back to the top.
% Type / followed by a keyword (like /error or /crash) and hit Enter to search. Press n to jump to the next match
% Press q to quit.

clear; close all; clc;
%% Cleanly reset the parallel environment before starting
% Order matters: kill workers first (so MathWorksServiceHost state is clean),
% THEN delete job files (rm-rf without killing first causes Index OOB on next parpool).
delete(gcp('nocreate'));
% my_pid = num2str(feature('getpid'));
% [~,~] = system(['pgrep -f "MATLAB/R2024b" | grep -v ^' my_pid '$ | xargs -r kill -9 2>/dev/null']);
% pause(1); % let workers fully die before touching job files
% matlab_dir = fileparts(prefdir);
% cluster_dir = fullfile(matlab_dir, 'local_cluster_jobs', ['R' version('-release')]);
% [~, ~] = system(['rm -rf ' cluster_dir '/Job* 2>/dev/null']);

%% Master Configuration
heights_km                          = 500:10:1200;
Master_config.Lat_range_deg         = [54+(35/60), 83+(40/60)]; %54°35N Denmark minimum, 83°40N Greenland max
% Master_config.Lat_range_deg         = [0 83.6];
Master_config.Min_elevation_UE      = 20;
Master_config.Num_Planes            = 2:25; % Num Planes
Master_config.Sats_Plane            = 2:25; % Sats per Plane
Master_config.Inc_vec               = linspace(70, 80, 111); % More general -> linspace(max(Lat_range_deg)-15, min(max(Lat_range_deg),80), 21)
Master_config.Target_num_candidates = 1;
Master_config.SampleTime            = 660;

% Sub Run configurations
Master_config.Ultrafast.Duration_h  = 3;  
Master_config.Ultrafast.Num_UEs     = 200;
Master_config.Fast.Duration_h       = 50;  
Master_config.Fast.Num_UEs          = 800;
certainty = 99 * 1e-2;
fractional_area = 0.1 * 1e-2;
fractional_time = 0.1 * 1e-2;
required_samples = log(1-certainty)/log(1-fractional_area*fractional_time)

required_time_h = ceil(sqrt(required_samples)) / (3600/Master_config.SampleTime)
Master_config.Detailed.Duration_h   = required_time_h;
required_UEs = ceil(sqrt(required_samples))
Master_config.Detailed.Num_UEs      = required_UEs;

Master_config.Max_workers = 16;

% Set resume_dir to a previous Master_Sweep folder to continue from where it
% left off, e.g. resume_dir = 'simulation_output/Master_Sweep_20260522_105447';
% Leave empty to start a fresh run.
resume_dir = 'simulation_output/Master_Sweep_20260530_135908';

% Record start time for the sweep
start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
fprintf('=======================================================\n');
fprintf('STARTING CONSTELLATION SWEEP AT: %s\n', char(start_time));
fprintf('=======================================================\n');

% Create master sweep output folder (or reuse existing one when resuming)
if ~isempty(resume_dir) && exist(resume_dir, 'dir')
    out_dir = resume_dir;
    checkpoint_file = fullfile(out_dir, 'Master_Altitude_Sweep_Results.mat');
    if ~exist(checkpoint_file, 'file')
        error('resume_dir specified but no checkpoint .mat found in: %s', out_dir);
    end
    ck = load(checkpoint_file);
    best_delta_sats = ck.best_delta_sats;
    star_sats       = ck.star_sats;

    % Build altitude -> gridsearch_dir map from the checkpoint.
    % Using a map avoids index-mismatch crashes when the sweep range changes
    % between runs (gridsearch_dirs and heights_km can have different lengths).
    n_ck = min(numel(ck.heights_km), numel(ck.gridsearch_dirs));
    ck_dir_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for ck_i = 1:n_ck
        if ~isempty(ck.gridsearch_dirs{ck_i})
            ck_dir_map(ck.heights_km(ck_i)) = ck.gridsearch_dirs{ck_i};
        end
    end

    % Re-build gridsearch_dirs and all_delta_sats indexed for the CURRENT heights_km
    heights_km_sorted = sort(heights_km, 'descend');
    gridsearch_dirs = cell(length(heights_km_sorted), 1);
    all_delta_sats  = cell(length(heights_km_sorted), 1);
    for ck_i = 1:length(heights_km_sorted)
        h = heights_km_sorted(ck_i);
        if isKey(ck_dir_map, h)
            gridsearch_dirs{ck_i} = ck_dir_map(h);
        end
    end

    completed_heights = heights_km_sorted(~cellfun(@isempty, gridsearch_dirs));
    fprintf('Resuming from: %s\n', out_dir);
    fprintf('Already completed %d altitudes. Skipping: %s km\n', ...
        numel(completed_heights), num2str(completed_heights));
else
    date_str_start = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    folder_name    = sprintf('Master_Sweep_%s', date_str_start);
    out_dir        = fullfile('simulation_output', folder_name);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    star_sats       = [];
    best_delta_sats = table();
    all_delta_sats  = {};
    gridsearch_dirs = cell(length(heights_km), 1);
    completed_heights = [];
end
gridsearch_base = fullfile(out_dir, 'gridsearch_runs');

% Runs through all orbit heights from top to bottom
heights_km = sort(heights_km, 'descend');
for i = 1:length(heights_km) 
    current_h_meters = heights_km(i) * 1000;
  
    %% Walker-Star calculation
    minimum_lat_deg = min(Master_config.Lat_range_deg);
    [Num_planes, Sats_per_plane, Total_sats] = calculate_walker_star(heights_km(i), minimum_lat_deg, Master_config.Min_elevation_UE);
    star_sats(i).Orbit_height   = heights_km(i);
    star_sats(i).Total_sats     = Total_sats;
    star_sats(i).Num_planes     = Num_planes;
    star_sats(i).Phasing        = Num_planes/2;
    star_sats(i).Inclination    = 90;
    star_sats(i).Sats_per_plane = Sats_per_plane;


    % Skip altitudes already completed in a previous (or current) run
    if ismember(heights_km(i), completed_heights)
        fprintf('Skipping already-completed altitude %d km\n', heights_km(i));
        continue;
    end

    if i > 1 && height(best_delta_sats) > 0
        min_sats = best_delta_sats.Total_sats(end); % Limit search space based on previous result
    else
        min_sats = 0;
    end

    [best_params, all_delta_sats{i}, gridsearch_dirs{i}] = gridsearch(Master_config, heights_km(i), min_sats, gridsearch_base);

    % Keep one best row per altitude in sweep order.
    best_delta_sats = [best_delta_sats; best_params(1, :)];

    % --- CHECKPOINT: save after every completed altitude ---
    save(fullfile(out_dir, 'Master_Altitude_Sweep_Results.mat'), ...
        'heights_km', 'star_sats', 'best_delta_sats', 'all_delta_sats', 'Master_config', 'gridsearch_dirs');
    fprintf('[Checkpoint] Saved after altitude %d km\n', heights_km(i));

    % Reset parallel environment between iterations: kill workers first, then purge files.
    % delete(gcp('nocreate'));
    % [~,~] = system(['pgrep -f "MATLAB/R2024b" | grep -v ^' my_pid '$ | xargs -r kill -9 2>/dev/null']);
    % pause(1);
    % [~, ~] = system(['rm -rf ' cluster_dir '/Job* 2>/dev/null']);

end

%% --- 3. Save Outputs ---
% (out_dir was already created before the loop)

% Save results for plot regeneration via plotting_scripts/plot_sweep.m
save(fullfile(out_dir, 'Master_Altitude_Sweep_Results.mat'), ...
    'heights_km', 'star_sats', 'best_delta_sats', 'all_delta_sats', 'Master_config', 'gridsearch_dirs');
fprintf('Results saved. Regenerate plots with: plot_sweep(''%s'')\n', out_dir);

% Save optimal constellations to workspace root so get_cfg.m can load them
save('optimal_constellations.mat', 'heights_km', 'star_sats', 'best_delta_sats');
fprintf('Optimal constellations saved to optimal_constellations.mat\n');

% --- RECORD END TIME ---
end_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
elapsed_time = end_time - start_time;

fprintf('\n=======================================================\n');
fprintf('SWEEP FINISHED AT: %s\n', char(end_time));
fprintf('TOTAL ELAPSED TIME: %s\n', char(elapsed_time));
fprintf('Master plot saved as Star_vs_Delta_Comparison.png\n');
fprintf('=======================================================\n');
