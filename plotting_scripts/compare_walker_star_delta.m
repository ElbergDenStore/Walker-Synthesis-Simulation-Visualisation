function compare_walker_star_delta()
clear; close all; clc;
% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "compare_walker_star_delta"
%% ===== SHARED SCENARIO SETTINGS =====
height_km        = 1000;
min_elevation_UE = 20;

Lat_range_deg = [54+(35/60), 83+(40/60)];
Lon_range_deg = [-(73+(10/60)), 32+(30/60)];

StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
StopTime  = datetime('3-Jun-2025 11:59:59', 'TimeZone', 'UTC'); % 48 hours

% Ku-band link budget
f_DL           = 12e9;
B_DL           = 50e6;
NF_DL          = 5;
G_rx           = 33;        % dBi, UE receive antenna gain
Target_PFD_MHz = -115;      % dBW/m²/MHz

FRF = 3;
RU  = 1;

NumUEs = 2000;

lat_vec = linspace(Lat_range_deg(1), Lat_range_deg(2), 8);       % 6 latitudes
lon_vec = linspace(Lon_range_deg(1), Lon_range_deg(2), 4);       % 6 longitudes
[UE_lats,UE_lons ] = meshgrid(lat_vec, lon_vec);

% [UE_lats, UE_lons] = generate_equal_area_ues(Lat_range_deg, Lon_range_deg, NumUEs);

%% ===== BUILD BASE CFG (fields shared by both constellations) =====
BaseCfg.Orbit_height             = height_km * 1e3;  % m
BaseCfg.Min_elevation_UE         = min_elevation_UE;
BaseCfg.SampleTime               = 60;               % s
BaseCfg.FRF                      = FRF;
BaseCfg.RU                       = RU;
BaseCfg.Use_P618                 = true;
BaseCfg.Modified_shannon         = true;
% BaseCfg.Simple_Atmospheric_Loss_dB = 1;
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

[Num_planes_star, Sats_per_plane_star, ~] = calculate_walker_star(height_km, min(Lat_range_deg), min_elevation_UE);
CfgStar.Num_planes     = Num_planes_star;
CfgStar.Sats_per_plane = Sats_per_plane_star;
CfgStar.Total_sats     = Num_planes_star * Sats_per_plane_star;
CfgStar.Inclination    = 90;            % Walker-Star: polar/near-polar
CfgStar.Phasing        = floor(Num_planes_star / 2);

fprintf('Walker Star: %d planes x %d sats/plane = %d total\n', ...
    CfgStar.Num_planes, CfgStar.Sats_per_plane, CfgStar.Total_sats);

%% ===== WALKER DELTA CONFIGURATION =====
CfgDelta = BaseCfg;
CfgDelta.WalkerStar = false;

CfgDelta.Num_planes     = 4;
CfgDelta.Sats_per_plane = 14;
CfgDelta.Inclination    = 77;
CfgDelta.Phasing        = 2;
CfgDelta.Total_sats     = CfgDelta.Num_planes * CfgDelta.Sats_per_plane;

fprintf('Walker Delta: %d planes x %d sats/plane = %d total\n', ...
    CfgDelta.Num_planes, CfgDelta.Sats_per_plane, CfgDelta.Total_sats);

%% ===== RUN SIMULATIONS =====
use_parallel = false;
calc_link    = true;

fprintf('\nRunning Walker Star simulation...\n');
metrics_star = constellation_simulator(CfgStar, use_parallel, calc_link);

fprintf('\nRunning Walker Delta simulation...\n');
metrics_delta = constellation_simulator(CfgDelta, use_parallel, calc_link);

%% ===== PLOT RESULTS =====
plot_simulation(metrics_star,  use_parallel);
plot_simulation(metrics_delta, use_parallel);

%% ===== SHOW CONSTELLATIONS =====
show_interactive = false;
save_fig = true;

if ~isempty(getenv('DISPLAY'))
    show_constellation(CfgStar,  show_interactive, save_fig);
    show_constellation(CfgDelta, show_interactive, save_fig);
else
    disp('Skipping show_constellation: no display available');
end

end
