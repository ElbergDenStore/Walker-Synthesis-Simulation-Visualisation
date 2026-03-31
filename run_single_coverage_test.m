clear all; close all; clc;

Cfg = get_cfg(1000,"walkerdelta","medium","medium");
Cfg.FRF = 3;
Cfg.RU = 1;

calc_link = true;
use_parallel = false;
metrics = coverage_simulator_function(Cfg,use_parallel,calc_link);

plot_simulation(metrics, use_parallel);

% save('last_run.mat', 'metrics')
% load('last_run.mat') 

% show_interactive = true;
% save_fig = false;
% 
% show_constellation(Cfg, show_interactive, save_fig)
% fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);
