% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run_single_coverage_test"
% /opt/VirtualGL/bin/vglrun matlab & - run it interactively with GPU rendering
function run_single_coverage_test()
clear all; close all; clc;

Cfg = get_cfg(1000,"walkerdelta","medium","medium");
Cfg.FRF = 3;
Cfg.RU = 1;

Cfg.DL.Direction     = "DL";
Cfg.DL.G_rx          = 39;  
Cfg.Target_PFD_MHz   = -128;
Cfg.DL.G_tx          = get_adjusted_tx_gain(Cfg.Orbit_height, Cfg.Min_elevation_UE, Cfg.DL.f); 
Cfg.DL.Max_P_tx_dBm  = PFD_calc(Cfg.Target_PFD_MHz, Cfg.DL.G_tx, Cfg.DL.B, Cfg.Orbit_height, Cfg.Min_elevation_UE);
Cfg.DL.Max_EIRP_dBm  = Cfg.DL.Max_P_tx_dBm + Cfg.DL.G_tx;

[Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons] = generate_equal_ish_area_UEs([54+(35/60), 83+(40/60)], [-(73+(10/60)), 33+(30/60)], 200);




calc_link = true;
use_parallel = true;
% metrics = fast_coverage_simulator_function(Cfg,false,false,false);

metrics = coverage_simulator_function(Cfg,use_parallel,calc_link);

Num_UEs = length(metrics.UEs);
plot_simulation(metrics, use_parallel);

% save('last_run.mat', 'metrics')
% load('last_run.mat') 

show_interactive = false;
save_fig = true;
% 
show_constellation(Cfg, show_interactive, save_fig)
% fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);
end