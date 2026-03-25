clear all; close all; clc;

Cfg = get_cfg(1000,"walkerdelta","big","medium");


calc_link = false;
plot_results = false;
if length(Cfg.Lat_vec)*length(Cfg.Lon_vec) > 8 % automate le decision
    use_parallel = true;
else
    use_parallel = false;
end
use_parallel = false;

metrics = coverage_simulator_function(Cfg, plot_results,use_parallel,calc_link);
show_interactive = true;
save_fig = false;

show_constellation(Cfg, show_interactive, save_fig)
fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);
