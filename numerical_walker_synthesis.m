% NUMERICAL_WALKER_SYNTHESIS  Press Run to search for minimal valid constellations.
%
%   This is the main entry point for the *numerical* (Monte-Carlo gridsearch)
%   constellation synthesis. For each orbital altitude it searches the
%   plane/satellite/inclination grid for the smallest constellation that
%   still meets the coverage requirement, then saves a resumable checkpoint.
%
%   HOW TO RUN
%     - Press the green Run button (or F5), or
%     - Headless overnight:
%         matlab -nosplash -nodesktop -batch "numerical_walker_synthesis" > log.txt
%       Follow along with:  tail -f log.txt
%
%   Edit the USER CONFIG block below, then run. Everything else is handled by
%   synthesis_preset() (defaults) and run_altitude_sweep() (the sweep engine).
% =========================================================================

clear; close all; clc;
delete(gcp('nocreate'));   % clean parallel state before a long run

% Ensure project root + all functions/ subfolders are on the path
repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'functions'));
path_setup();

%% ===================== USER CONFIG ======================================
% 1) Pick a preset. It fills Master_config + a default altitude vector.
%      "regional_delta" - Walker Delta, Arctic/Nordic (54.6-83.7 N), i swept 70-80
%      "global_delta"   - Walker Delta, tropical band (0-30 N), wider grid
%      "star"           - Walker Star, polar (i = 90), Arctic/Nordic
preset = "regional_delta";

% 2) Resume a previous run? Point at its folder to continue where it stopped.
%    Leave '' to start fresh.
%      e.g. resume_dir = 'simulation_output/Master_Sweep_20260530_135908';
resume_dir = '';

[Master_config, heights_km] = synthesis_preset(preset);

% 3) Common knobs to override (delete/comment any line to keep the preset default).
%    heights_km: a vector sweeps; a scalar runs a single altitude.
% heights_km                   = 500:10:1200;          % [km]
% Master_config.Lat_range_deg  = [54+(35/60), 83+(40/60)];  % [min max] latitude band
% Master_config.Min_elevation_UE = 20;                 % [deg] min elevation
% Master_config.Num_Planes     = 2:25;                 % planes to search
% Master_config.Sats_Plane     = 2:25;                 % sats/plane to search
% Master_config.Inc_vec        = linspace(70, 80, 111);% [deg] inclinations to search
% Master_config.Max_workers    = 16;                   % parallel pool size
%% =======================================================================

out_dir = run_altitude_sweep(Master_config, heights_km, resume_dir);
fprintf('\nDone. Sweep output: %s\n', out_dir);
