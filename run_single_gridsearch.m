% matlab -nosplash -nodesktop -batch "run_single_gridsearch" > "log.txt"
clear; close all; clc;
% Kill workers first, THEN delete job files.
% Deleting files while workers are alive leaves MathWorksServiceHost with stale
% in-memory references, causing 'Index exceeds array bounds' on the next parpool().
delete(gcp('nocreate'));
my_pid = num2str(feature('getpid'));
[~,~] = system(['pgrep -f "MATLAB/R2024b" | grep -v ^' my_pid '$ | xargs -r kill -9 2>/dev/null']);
pause(1); % let workers fully die
matlab_dir = fileparts(prefdir);
cluster_dir = fullfile(matlab_dir, 'local_cluster_jobs', ['R' version('-release')]);
[~, ~] = system(['rm -rf ' cluster_dir '/Job* 2>/dev/null']);


Master_config.Lat_range_deg         = [54+(35/60), 83+(40/60)]; %54°35N Denmark minimum, 83°40N Greenland max
Master_config.Min_elevation_UE      = 20;
Master_config.Num_Planes            = 2:32; % Num Planes
Master_config.Sats_Plane            = 2:32; % Sats per Plane
Master_config.Inc_vec               = linspace(70, 80, 81); % More general -> linspace(max(Lat_range_deg)-15, min(max(Lat_range_deg),80), 21)
Master_config.Target_num_candidates = 1;
Master_config.SampleTime            = 420;

% Sub Run configurations
Master_config.Ultrafast.Duration_h  = 2;  
Master_config.Ultrafast.Num_UEs     = 200;
Master_config.Fast.Duration_h       = 40;  
Master_config.Fast.Num_UEs          = 800;
certainty = 99 * 1e-2;
fractional_area = 0.1 * 1e-2;
fractional_time = 0.1 * 1e-2;
required_samples = log(1-certainty)/log(1-fractional_area*fractional_time)

required_time_h = ceil(sqrt(required_samples)) / (3600/Master_config.SampleTime)
Master_config.Detailed.Duration_h   = required_time_h;
required_UEs = ceil(sqrt(required_samples))
Master_config.Detailed.Num_UEs      = required_UEs;



orbit_height_km = 1000;
min_sats = 55;
plot_individual_results = true;
[best_params, all_delta_sats] = gridsearch(Master_config, orbit_height_km, plot_individual_results, min_sats);
