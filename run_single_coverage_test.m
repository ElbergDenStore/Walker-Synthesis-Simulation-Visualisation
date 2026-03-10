clear all; close all; clc;

Cfg.StartTime  = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
Cfg.StopTime   = datetime('1-Jun-2025 0:59:59', 'TimeZone', 'UTC');
Cfg.SampleTime = 20; % seconds

Cfg.Lat_vec = linspace(55, 85, 2); 
Cfg.Lon_vec = linspace(-60, 30, 2);

% Constellation
Cfg.Min_elevation_UE = 20;
% Cfg.Orbit_height = 932e3;
% Cfg.Num_planes     = 5;  
% Cfg.Sats_per_plane = 14;
Cfg.Orbit_height = 707e3;
Cfg.Num_planes     = 6;  
Cfg.Sats_per_plane = 17;
Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
Cfg.Inclination    = 88; 
Cfg.Phasing = Cfg.Num_planes / 2;
Cfg.WalkerStar     = true;  
Normal_gap         = 180 / (Cfg.Num_planes - 1/3); % only valid for flattop hexagons
Cfg.Seam_gap       = 2/3 * Normal_gap;             % only valid for flattop hexagons


% Downlink Link Budget Config FR2
Cfg.DL.Direction = "DL";
Cfg.DL.B         = 2e6;     
Cfg.DL.f         = 20e9;    
Cfg.DL.P_tx_dBm  = 20; 
Cfg.DL.G_tx      = 42; 
Cfg.DL.Tx_type   = "array";      
Cfg.DL.G_rx      = 32;      
Cfg.DL.Rx_type   = "array";      
Cfg.DL.NF        = 5;
Cfg.DL.EIRP_dBm  = Cfg.DL.P_tx_dBm + Cfg.DL.G_tx;

% % Downlink Link Budget Config FR1
% Cfg.DL.Direction = "DL";
% Cfg.DL.B         = 180e3;
% Cfg.DL.f         = 2.6e9;
% Cfg.DL.P_tx_dBm  = 10 + 30 + 10*log10(Cfg.DL.B/1e6); %EIRP 34 dBW /MHz TR 38821 %26.6 dBm @ 180kHz
% Cfg.DL.G_tx      = 24; %TR 38821
% Cfg.DL.Tx_type   = "array";      
% Cfg.DL.G_rx      = 3; % Dipole antenna typically
% Cfg.DL.Rx_type   = "array";      
% Cfg.DL.NF        = 7; %  TR 38821
% Cfg.DL.EIRP_dBm  = Cfg.DL.P_tx_dBm + Cfg.DL.G_tx;



calc_link = false;
plot_results = false;
if length(Cfg.Lat_vec)*length(Cfg.Lon_vec) > 8 % automate le decision
    use_parallel = true;
else
    use_parallel = false;
end

metrics = coverage_simulator_function(Cfg, plot_results,use_parallel,calc_link);
fprintf("Worst Coverage percentage" + metrics.worst_coverage_percent);