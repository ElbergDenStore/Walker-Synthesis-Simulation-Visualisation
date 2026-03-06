clear all; close all; clc;

Cfg.StartTime  = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
Cfg.StopTime   = datetime('1-Jun-2025 2:59:59', 'TimeZone', 'UTC');
Cfg.SampleTime = 20; % seconds

Cfg.Lat_vec = linspace(50, 85, 2); 
Cfg.Lon_vec = linspace(-60, 30, 2);

% Constellation
Cfg.Min_elevation_UE = 20;
Cfg.Orbit_height = 932e3;
Cfg.Num_planes     = 5;  
Cfg.Sats_per_plane = 14;
Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
Cfg.Inclination    = 88; 
Cfg.Phasing = Cfg.Num_planes / 2;
Cfg.WalkerStar     = true;  
Normal_gap         = 180 / (Cfg.Num_planes - 1/3); % only valid for flattop hexagons
Cfg.Seam_gap       = 2/3 * Normal_gap;             % only valid for flattop hexagons


% Downlink Link Budget Config
Cfg.DL.Direction = "DL";
Cfg.DL.B         = 2e6;     
Cfg.DL.f         = 26e9;    
Cfg.DL.P_tx_dBm  = 30; 
Cfg.DL.G_tx      = 36; 
Cfg.DL.Tx_type   = "array";      
Cfg.DL.G_rx      = 38.3;      
Cfg.DL.Rx_type   = "dish";      
Cfg.DL.NF        = 5;
Cfg.DL.EIRP_dBm  = Cfg.DL.P_tx_dBm + Cfg.DL.G_tx;


plot_results = true;
if length(Cfg.Lat_vec)*length(Cfg.Lon_vec) > 8 % automate le decision
    use_parallel = true;
else
    use_parallel = false;
end

metrics = coverage_simulator_function(Cfg, plot_results,use_parallel);
fprintf("Worst Coverage percentage" + metrics.worst_coverage_percent);