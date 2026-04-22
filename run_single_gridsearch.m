clear; close all; clc;
delete(gcp('nocreate')); % necessary or it might get stuck


Master_config.Lat_range_deg         = [55+(35/60), 83+(40/60)]; %54°35N Denmark minimum, 83°40N Greenland max
Master_config.Min_elevation_UE      = 20;
Master_config.Num_Planes            = 2:20; % Num Planes
Master_config.Sats_Plane            = 2:20; % Sats per Plane
Master_config.Inc_vec               = linspace(70, 80, 41); % More general -> linspace(max(Lat_range_deg)-15, min(max(Lat_range_deg),80), 21)
Master_config.Target_num_candidates = 1000;

% Sub Run configurations
Master_config.Ultrafast.Duration_h  = 1;  
Master_config.Ultrafast.Num_UEs     = 100;
Master_config.Fast.Duration_h       = 12;  
Master_config.Fast.Num_UEs          = 300;
Master_config.Detailed.Duration_h   = 36;      
Master_config.Detailed.Num_UEs      = 2000;


orbit_height_km = 1000;
min_sats = 0;
plot_individual_results = true;
[best_params, all_delta_sats{i}] = gridsearch(Master_config, orbit_height_km, plot_individual_results, min_sats);
