clear all; close all; clc;

Cfg = get_cfg(1000,"walkerdelta","big","medium");
Cfg.FRF = 3;
Cfg.RU = 1;

calc_link = true;
plot_results = false;
use_parallel = false;
metrics = coverage_simulator_function(Cfg, plot_results,use_parallel,calc_link);


% show_interactive = true;
% save_fig = false;
% 
% show_constellation(Cfg, show_interactive, save_fig)
% fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);
