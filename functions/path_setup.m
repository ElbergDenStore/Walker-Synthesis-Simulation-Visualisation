function path_setup()
% PATH_SETUP  Add the project root and ALL functions/ subfolders to the path.
%   This file lives in functions/, so the project root is one folder up.
%   genpath() is used so any nested subfolder of functions/ (e.g.
%   functions/link_budget/, functions/generators/) is resolved automatically.
%   You can move helper files into functions/<subfolder>/ freely without
%   breaking MATLAB's filename-based function lookup.
    this_dir     = fileparts(mfilename('fullpath'));   % .../functions
    project_root = fileparts(this_dir);                % repo root
    addpath(project_root);
    addpath(genpath(fullfile(project_root, 'functions')));  % recursive; covers functions/data
end