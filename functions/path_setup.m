function path_setup()
    % Find where this setup file lives (inside plotting_scripts)
    [this_dir, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Define the project root as one folder up
    project_root = fullfile(this_dir, '..');
    
    % Add the paths to the MATLAB session
    addpath(project_root);
    addpath(fullfile(project_root, 'functions'));
    addpath(fullfile(project_root, 'functions', 'data'));
end