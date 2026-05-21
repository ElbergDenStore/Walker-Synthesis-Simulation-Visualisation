% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run_single_coverage_test"
% /opt/VirtualGL/bin/vglrun matlab & - run it interactively with GPU rendering
function run_single_coverage_test()
clear all; close all; clc;
% Cfg = get_cfg(1000,"walkerdelta","medium","long");

height_km        = 510;
min_elevation_UE = 20;

Lat_range_deg = [54+(35/60), 83+(40/60)];
% Lon_range_deg = [-(73+(10/60)), 33+(30/60)];
Lon_range_deg = [-180, 180];
StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
StopTime  = datetime('1-Jun-2025 12:59:59', 'TimeZone', 'UTC'); % 48 hours

% Ku-band link budget
f_DL           = 12e9;
B_DL           = 50e6;
NF_DL          = 5;
G_rx           = 33;        % dBi, UE receive antenna gain
Target_PFD_MHz = -115;      % dBW/m²/MHz

FRF = 3;
RU  = 1;

NumUEs = 200;

[UE_lats, UE_lons] = generate_equal_ish_area_UEs(Lat_range_deg, Lon_range_deg, NumUEs);

%% ===== BUILD BASE CFG (fields shared by both constellations) =====
BaseCfg.Orbit_height             = height_km * 1e3;  % m
BaseCfg.Min_elevation_UE         = min_elevation_UE;
BaseCfg.SampleTime               = 60;               % s
BaseCfg.FRF                      = FRF;
BaseCfg.RU                       = RU;
BaseCfg.Use_P618                 = false;
BaseCfg.Modified_shannon         = true;
BaseCfg.Simple_Atmospheric_Loss_dB = 1;
BaseCfg.Share_bandwidth          = false;
BaseCfg.Target_PFD_MHz           = Target_PFD_MHz;
BaseCfg.StartTime                = StartTime;
BaseCfg.StopTime                 = StopTime;
BaseCfg.Flat_UE_array.Lats       = UE_lats;
BaseCfg.Flat_UE_array.Lons       = UE_lons;

BaseCfg.DL.Direction     = "DL";
BaseCfg.DL.f             = f_DL;
BaseCfg.DL.B             = B_DL;
BaseCfg.DL.NF            = NF_DL;
BaseCfg.DL.G_rx          = G_rx;
BaseCfg.DL.Tx_type       = "array";
BaseCfg.DL.Rx_type       = "array";
BaseCfg.DL.G_tx          = get_adjusted_tx_gain(BaseCfg.Orbit_height, min_elevation_UE, f_DL);
BaseCfg.DL.Max_P_tx_dBm  = PFD_calc(Target_PFD_MHz, BaseCfg.DL.G_tx, B_DL, BaseCfg.Orbit_height, min_elevation_UE);
BaseCfg.DL.Max_EIRP_dBm  = BaseCfg.DL.Max_P_tx_dBm + BaseCfg.DL.G_tx;
BaseCfg.DL.Max_EIRP_dBm_Hz = BaseCfg.DL.Max_EIRP_dBm - 10*log10(B_DL);
BaseCfg.DL.BeamGrid      = calculate_hexagonal_beams(BaseCfg.DL.G_tx, f_DL, BaseCfg.Orbit_height, min_elevation_UE, BaseCfg.DL.Max_EIRP_dBm_Hz, FRF);

%% ===== WALKER STAR CONFIGURATION =====
CfgStar = BaseCfg;
CfgStar.WalkerStar = true;

% 4               14              90            2           56             1100      
% [Num_planes_star, Sats_per_plane_star, ~] = calculate_walker_star(height_km, min(Lat_range_deg), min_elevation_UE);
CfgStar.Num_planes     = 7;
CfgStar.Sats_per_plane = 23;
CfgStar.Total_sats     = CfgStar.Num_planes * CfgStar.Sats_per_plane;
CfgStar.Inclination    = 90; 
CfgStar.Phasing        = CfgStar.Num_planes / 2;

%% ===== RUN SIMULATIONS =====
use_parallel = false;
calc_link    = false;

% fprintf('\nRunning Walker Star simulation...\n');
metrics_star = coverage_simulator_function(CfgStar, use_parallel, calc_link);
fprintf("coverage percentage %0.8f",metrics_star.worst_coverage_percent)

%% ===== PLOT RESULTS =====
% plot_simulation(metrics_star,  use_parallel);

show_interactive = true;
save_fig = false;
show_constellation(CfgStar, show_interactive, save_fig)

end