% COMPARE_UE_DISTRIBUTION  Compare uniform vs population-based UE distributions.
%   Runs the same Walker Star constellation twice - once with equal-area UEs and
%   once with population-weighted UEs - over the same scenario, and plots both
%   coverage results. Edit the shared parameters below, then run.
%
% Headless run:
%   xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nodesktop -batch "compare_ue_distribution"
clear all; close all; clc;
delete(gcp('nocreate'));
% Add the project to the MATLAB path (robust to the script's folder depth).
repo_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
addpath(fullfile(repo_root, 'functions'));
path_setup();
%% ===== SHARED PARAMETERS =====
height_km        = 1000;
min_elevation_UE = 20;

Lat_range_deg = [54+(35/60), 83+(40/60)];
Lon_range_deg = [-(73+(10/60)), 33+(30/60)];
StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
StopTime  = datetime('1-Jun-2025 14:59:59', 'TimeZone', 'UTC');

% Ku-band link budget
f_DL           = 12e9;
B_DL           = 250e6;
NF_DL          = 5;
G_rx           = 33;        % dBi, UE receive antenna gain
Target_PFD_MHz = -123;      % dBW/m²/MHz

FRF = 3;
RU  = 1;

% UE generation parameters
people_per_ue    = 3000;   % for population-based distribution (uniform count is derived to match)

%% ===== WALKER STAR CONSTELLATION =====
Num_planes     = 5;
Sats_per_plane = 13;
Inclination    = 90;

%% ===== BUNDLE SHARED PARAMETERS FOR build_cfg =====
P = struct('height_km', height_km, 'min_elevation_UE', min_elevation_UE, ...
           'FRF', FRF, 'RU', RU, 'Target_PFD_MHz', Target_PFD_MHz, ...
           'StartTime', StartTime, 'StopTime', StopTime, ...
           'f_DL', f_DL, 'B_DL', B_DL, 'NF_DL', NF_DL, 'G_rx', G_rx, ...
           'Num_planes', Num_planes, 'Sats_per_plane', Sats_per_plane, ...
           'Inclination', Inclination);

use_parallel = false;
calc_link    = true;

%% ===== SHARED OUTPUT DIRECTORY =====
script_dir = fileparts(mfilename('fullpath'));
date_str   = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
run_dir    = fullfile(script_dir, '..', 'figures', 'ue_distribution', sprintf('Compare_UE_%s', date_str));
mkdir(fullfile(run_dir, 'uniform'));
mkdir(fullfile(run_dir, 'population'));
fprintf('Saving all outputs to: %s\n', run_dir);

%% ===== GENERATE POPULATION-BASED UEs FIRST (to determine matching UE count) =====
fprintf('Generating population-based UEs to determine UE count...\n');
[UE_lats_p, UE_lons_p, Total_Pop] = generate_population_ues(Lat_range_deg, Lon_range_deg, people_per_ue);
NumUEs_uniform = length(UE_lats_p);
fprintf('UE count for this region: %d (will be used for both distributions)\n', NumUEs_uniform);

%% ===== RUN: UNIFORM (EQUAL-AREA) DISTRIBUTION =====
fprintf('\n========== [1/2] Uniform UE Distribution ==========\n');
[UE_lats_u, UE_lons_u] = generate_equal_area_ues(Lat_range_deg, Lon_range_deg, NumUEs_uniform);
CfgUniform = build_cfg(P, UE_lats_u, UE_lons_u);
CfgUniform.Total_Pop = 0;  % no population data for uniform
CfgUniform.Save_dir  = fullfile(run_dir, 'uniform');

metrics_uniform = Constellation_simulator(CfgUniform, use_parallel, calc_link);
fprintf('Uniform   -> worst coverage: %.4f%%\n', metrics_uniform.worst_coverage_percent);
plot_simulation(metrics_uniform, use_parallel);

%% ===== RUN: POPULATION-BASED DISTRIBUTION =====
fprintf('\n========== [2/2] Population-Based UE Distribution ==========\n');
CfgPop = build_cfg(P, UE_lats_p, UE_lons_p);
CfgPop.Total_Pop = Total_Pop;
CfgPop.Save_dir  = fullfile(run_dir, 'population');

metrics_pop = Constellation_simulator(CfgPop, use_parallel, calc_link);
fprintf('Population -> worst coverage: %.4f%%\n', metrics_pop.worst_coverage_percent);
plot_simulation(metrics_pop, use_parallel);

fprintf('\nDone. Both simulation outputs saved to simulation_output/\n');

%% ===== BUILD COMMON LINK BUDGET FIELDS =====
function Cfg = build_cfg(P, UE_lats, UE_lons)
    Cfg.Orbit_height               = P.height_km * 1e3;
    Cfg.Min_elevation_UE           = P.min_elevation_UE;
    Cfg.SampleTime                 = 60;
    Cfg.FRF                        = P.FRF;
    Cfg.RU                         = P.RU;
    Cfg.Use_P618                   = false;
    Cfg.Modified_shannon           = true;
    Cfg.Simple_Atmospheric_Loss_dB = 1;
    Cfg.Share_bandwidth            = true;
    Cfg.Target_PFD_MHz             = P.Target_PFD_MHz;
    Cfg.StartTime                  = P.StartTime;
    Cfg.StopTime                   = P.StopTime;
    Cfg.Flat_UE_array.Lats         = UE_lats;
    Cfg.Flat_UE_array.Lons         = UE_lons;

    Cfg.DL.Direction    = "DL";
    Cfg.DL.f            = P.f_DL;
    Cfg.DL.B            = P.B_DL;
    Cfg.DL.NF           = P.NF_DL;
    Cfg.DL.G_rx         = P.G_rx;
    Cfg.DL.Tx_type      = "array";
    Cfg.DL.Rx_type      = "array";
    Cfg.DL.G_tx         = get_adjusted_tx_gain(Cfg.Orbit_height, P.min_elevation_UE, P.f_DL);
    Cfg.DL.Max_P_tx_dBm = PFD_calc(P.Target_PFD_MHz, Cfg.DL.G_tx, P.B_DL, Cfg.Orbit_height, P.min_elevation_UE);
    Cfg.DL.Max_EIRP_dBm = Cfg.DL.Max_P_tx_dBm + Cfg.DL.G_tx;
    Cfg.DL.Max_EIRP_dBm_Hz = Cfg.DL.Max_EIRP_dBm - 10*log10(P.B_DL);
    Cfg.DL.BeamGrid     = calculate_hexagonal_beams(Cfg.DL.G_tx, P.f_DL, Cfg.Orbit_height, P.min_elevation_UE, Cfg.DL.Max_EIRP_dBm_Hz, P.FRF);

    Cfg.WalkerStar      = true;
    Cfg.Num_planes      = P.Num_planes;
    Cfg.Sats_per_plane  = P.Sats_per_plane;
    Cfg.Total_sats      = P.Num_planes * P.Sats_per_plane;
    Cfg.Inclination     = P.Inclination;
    Cfg.Phasing         = P.Num_planes / 2;
end
