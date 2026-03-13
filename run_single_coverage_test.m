clear all; close all; clc;

Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
Cfg.StopTime   = datetime('2-Jun-2025 11:59:59', 'TimeZone', 'UTC');
Cfg.SampleTime = 60; % seconds

Cfg.Lat_vec = linspace(55, 85, 2); 
Cfg.Lon_vec = linspace(-60, 30, 2);

% Constellation
Cfg.Min_elevation_UE = 20;
% Cfg.Orbit_height = 932e3;
% Cfg.Num_planes     = 5;  
% Cfg.Sats_per_plane = 14;
Cfg.Orbit_height = 550e3;
Cfg.Num_planes     = 4;  
Cfg.Sats_per_plane = 14;
Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
Cfg.Inclination    = 77; 
% Cfg.Phasing = Cfg.Num_planes / 2;
Cfg.Phasing = 2;
Cfg.WalkerStar     = false;  
% Normal_gap         = 180 / (Cfg.Num_planes - 1/3); % only valid for flattop hexagons
% Cfg.Seam_gap       = 2/3 * Normal_gap;             % only valid for flattop hexagons


% Downlink Link Budget Config FR2
% PFD_regulation = -105;
Cfg.Target_PFD_MHz = -108;
Cfg.DL.Direction = "DL";
Cfg.DL.B         = 2e6;     
Cfg.DL.f         = 20e9;    
% Cfg.DL.P_tx_dBm  = 20; 
Cfg.DL.G_tx      = 42; 
Cfg.DL.Tx_type   = "array";      
Cfg.DL.G_rx      = 32;      
Cfg.DL.Rx_type   = "array";      
Cfg.DL.NF        = 5;

% attempt to change G based on orbit height to more directly compare the
% constellations
% Some reference
h_ref = 1200e3; %
G_tx_ref = 42;  % 
Re = 6371e3;    %

% Calculate Slant Range for the Reference Altitude
slant_ref = -Re*sind(Cfg.Min_elevation_UE) + sqrt(Re^2*sind(Cfg.Min_elevation_UE)^2 - (Re^2-(Re+h_ref)^2));

% Calculate Slant Range for the CURRENT Altitude in the sweep
current_h = Cfg.Orbit_height;
slant_current = -Re*sind(Cfg.Min_elevation_UE) + sqrt(Re^2*sind(Cfg.Min_elevation_UE)^2 - (Re^2-(Re+current_h)^2));

% Scale the Antenna Gain to keep the Ground Footprint constant
Cfg.DL.G_tx = G_tx_ref + 20 * log10(slant_current / slant_ref);


Cfg.DL.Max_P_tx_dBm  = PFD_calc(Cfg.Target_PFD_MHz, Cfg.DL.G_tx, Cfg.DL.B, Cfg.Orbit_height, Cfg.Min_elevation_UE);
Cfg.DL.Max_EIRP_dBm  = Cfg.DL.Max_P_tx_dBm + Cfg.DL.G_tx;

% % Downlink Link Budget Config FR1
% PFD_regulation = -113;
% Target_PFD = -116;
% Cfg.DL.Direction = "DL";
% Cfg.DL.B         = 180e3;
% Cfg.DL.f         = 2.6e9;
% Cfg.DL.G_tx      = 24; %TR 38821
% Cfg.DL.Tx_type   = "array";      
% Cfg.DL.G_rx      = 3; % Dipole antenna typically
% Cfg.DL.Rx_type   = "array";      
% Cfg.DL.NF        = 7; %  TR 38821
% Cfg.DL.P_tx_dBm  = 10;
% Cfg.DL.P_tx_dBm  = PFD_calc(Target_PFD_MHz, Cfg.DL.G_tx, Cfg.DL.B, Cfg.Orbit_height, Cfg.Min_elevation_UE);
% Cfg.DL.EIRP_dBm  = Cfg.DL.P_tx_dBm + Cfg.DL.G_tx;


calc_link = false;
plot_results = false;
if length(Cfg.Lat_vec)*length(Cfg.Lon_vec) > 8 % automate le decision
    use_parallel = true;
else
    use_parallel = false;
end

metrics = coverage_simulator_function(Cfg, plot_results,use_parallel,calc_link);
show_interactive = false;
save_fig = true;

show_constellation(Cfg, show_interactive, save_fig)
fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);
