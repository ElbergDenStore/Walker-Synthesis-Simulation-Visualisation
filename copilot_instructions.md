Software development guideline: 
Allow it to crash if the "protection code" is too complex. If I am finding all the valid indices by using >20 filter, it is possible to not have any valid indices, but only if the constellation is completely retarded or the simulation settings are fucked. It needs to be very human readable

Software style guide:
I need to be consistent, I dont care about following some IEEE standard or whatever, but i want it consistent.
Cfg.Num_sats
for structs and similar, the first letter needs to be capitalized. Cfg.Orbit_height_m instead of Cfg.orbit_height_m. if not part of a struct, orbit_height_m instead of Orbit_height_m. 
PFD is an acronym for Power Flux Density. I want functions and variables to have PFD in capital letters. dBm means dB and here i want the B to be capital. Same for Gbps and MHz and similar
I should have config instead of cfg to avoid confusion. Same for other things that should not be shortene? Num Sats versus Number_Satellites? Shit i dont know what is best...
If the variable unit can be misunderstood, it needs to have a post_fix or whatever orbit_height_m, gain_dBi p_dBm and so on.
The main simulators needs to be starting with capital letter instead of the current naming
I hate the naming of AI generated code such as "%% 2. The Smart Filters". It should just be "% Filters" If I see numbering I get angry.


I need the technical code to be clearly visible so noise needs to be moved to local functions: example:

script_dir     = fileparts(mfilename('fullpath'));
workspace_root = fileparts(script_dir);
addpath(fullfile(workspace_root, 'functions'));  % bootstrap so path_setup is found
path_setup();                                    % add repo root + all functions/ subfolders

sim_dir   = fullfile(workspace_root, 'simulation_output');
runs_root = fullfile(sim_dir, 'gridsearch_runs');

if nargin == 0
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

%%%% ACTUAL FUNCTIONALITY STARTS AFTER THIS %%%%

This is not nice, but the code is really needed for it to be easy to use.


Plots:
The plots made by the simulator needs to be consistent as well. 
Blue means low, yellow means high
Units always presented in parantheses "Orbit Height (m)" (Capitalize all important words)
The plots should be self contained. If the simulation is using beam utilization of 20% and that changes the results by 5x, it needs to be a part of the title. If it does not change anything, the parameter should not be shown.
The plots should be saved and report ready.

Comments in code:
Docstring in all files
I dont want any stupid comments like "-FIXED fast matrix computation now", but i do want to help people going through the source code explaining what is going on.




