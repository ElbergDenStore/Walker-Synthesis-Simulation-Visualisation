% matlab -nosplash -nodesktop -batch "run_single_gridsearch" > "log.txt"
clear; close all; clc;
delete(gcp('nocreate'));



Master_config.Lat_range_deg         = [54+(35/60), 83+(40/60)]; %54°35N Denmark minimum, 83°40N Greenland max
% Master_config.Lat_range_deg         = [34.5, 60.2]; % Europe from Cape Trypiti, Gavdos Island, Greece (34° 48′ 02″ N) to Helsinki 60°10′15″N 24°56′15″E
% Master_config.Lat_range_deg         = [0, 90];
Master_config.Min_elevation_UE      = 20;
Master_config.Num_Planes            = 2:25; % Num Planes
Master_config.Sats_Plane            = 2:35; % Sats per Plane
% Master_config.Inc_vec               = linspace(min(Master_config.Lat_range_deg), min(max(Master_config.Lat_range_deg),80), 81);
Master_config.Inc_vec               = linspace(70, 80, 111); % More general -> linspace(max(Lat_range_deg)-15, min(max(Lat_range_deg),80), 21)
Master_config.WalkerStar            = false;

% Master_config.Num_Planes            = 2:15;
% Master_config.Sats_Plane            = 5:35;
% Master_config.Inc_vec               = 90;
% Master_config.WalkerStar            = true;

Master_config.Target_num_candidates = 3;
Master_config.SampleTime            = 660;

% Sub Run configurations
Master_config.Ultrafast.Duration_h  = 3;  
Master_config.Ultrafast.Num_UEs     = 200;
Master_config.Fast.Duration_h       = 50;  
Master_config.Fast.Num_UEs          = 800;
certainty = 99 * 1e-2;
fractional_area = 0.1 * 1e-2;
fractional_time = 0.1 * 1e-2;
required_samples = log(1-certainty)/log(1-fractional_area*fractional_time)

required_time_h = ceil(sqrt(required_samples)) / (3600/Master_config.SampleTime)
Master_config.Detailed.Duration_h   = required_time_h;
required_UEs = ceil(sqrt(required_samples))
Master_config.Detailed.Num_UEs      = required_UEs;

Master_config.Max_workers = 2;

% orbit_height_km = 8000;
% min_sats = 18;
% Master_config.Min_elevation_UE      = 28;

orbit_height_km = 1000;
min_sats = 50;
% Master_config.Min_elevation_UE      = 30; %Global Delta
% Master_config.Min_elevation_UE      = 38; %Global Star
% Master_config.Min_elevation_UE      = 45; %Regional LEO Delta + Global MEO
% Master_config.Min_elevation_UE      = 43; %Regional LEO Star + Global MEO

[best_params, all_delta_sats] = gridsearch(Master_config, orbit_height_km, min_sats);
