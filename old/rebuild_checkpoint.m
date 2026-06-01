% rebuild_checkpoint.m  -- Reconstruct Master_Altitude_Sweep_Results.mat
% from whatever gridsearch_runs folders actually exist on disk.
% Usage: set sweep_dir below, then run.

sweep_dir = 'simulation_output/Master_Sweep_WalkerStar_20260522_095635';
runs_dir  = fullfile(sweep_dir, 'gridsearch_runs');

%% Scan all completed run folders
entries = dir(runs_dir);
entries = entries([entries.isdir] & ~startsWith({entries.name}, '.'));

% For each altitude keep only the LATEST folder (highest timestamp string)
alt_to_dir = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:numel(entries)
    tok = regexp(entries(k).name, '^(\d+)_\d+_(\d{8}_\d{6})$', 'tokens');
    if isempty(tok), continue; end
    alt = str2double(tok{1}{1});
    ts  = tok{1}{2};
    fp  = fullfile(runs_dir, entries(k).name);
    if ~isKey(alt_to_dir, alt)
        alt_to_dir(alt) = fp;
    else
        prev_tok = regexp(alt_to_dir(alt), '_(\d{8}_\d{6})$', 'tokens');
        if ~isempty(prev_tok) && str2double(strrep(ts,'_','')) > str2double(strrep(prev_tok{1}{1},'_',''))
            alt_to_dir(alt) = fp;
        end
    end
end

heights_km      = sort(cell2mat(keys(alt_to_dir)), 'descend');
gridsearch_dirs = cell(length(heights_km), 1);
all_star_sats   = cell(length(heights_km), 1);
best_star_sats  = table();
star_sats       = [];

for ii = 1:length(heights_km)
    h = heights_km(ii);
    gridsearch_dirs{ii} = alt_to_dir(h);

    pd_file = fullfile(alt_to_dir(h), 'plot_data.mat');
    if ~exist(pd_file, 'file'), continue; end
    pd = load(pd_file);
    cands = pd.plot_data.all_candidates;
    if isempty(cands), continue; end
    cands = sortrows(cands, 'Total_sats', 'ascend');
    best_row = cands(1, :);
    if ~ismember('Orbit_height_km', best_row.Properties.VariableNames)
        best_row.Orbit_height_km = h;
    end
    best_star_sats = [best_star_sats; best_row];
end

%% Load Master_config from existing mat (only field we can't reconstruct)
existing = load(fullfile(sweep_dir, 'Master_Altitude_Sweep_Results.mat'), 'Master_config');
Master_config = existing.Master_config;

%% Save
matfile_path = fullfile(sweep_dir, 'Master_Altitude_Sweep_Results.mat');
save(matfile_path, 'heights_km', 'star_sats', 'best_star_sats', ...
    'all_star_sats', 'Master_config', 'gridsearch_dirs');

fprintf('Checkpoint rebuilt with %d completed altitudes:\n', length(heights_km));
fprintf('  %d', heights_km); fprintf('\n');
disp(best_star_sats(:, {'Orbit_height_km','Total_sats','Num_planes','Sats_per_plane'}));

