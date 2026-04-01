clear all; close all; clc;

Cfg = get_cfg(1000,"walkerdelta","medium","medium");
Cfg.FRF = 3;
Cfg.RU = 1;

Cfg.Lat_vec = linspace(55, 85, 20); 
Cfg.Lon_vec = linspace(-180, 180, 3);

calc_link = true;
use_parallel = false;
metrics = coverage_simulator_function(Cfg,use_parallel,calc_link);
Num_UEs = length(metrics.UEs)
% plot_simulation(metrics, use_parallel);

% save('last_run.mat', 'metrics')
% load('last_run.mat') 

show_interactive = true;
save_fig = false;

show_constellation(Cfg, show_interactive, save_fig)
% fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);
