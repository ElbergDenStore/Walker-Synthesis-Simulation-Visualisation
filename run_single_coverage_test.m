clear all; close all; clc;

Cfg = get_cfg(1000,"walkerdelta","big","medium");

Cfg.Use_P618 = true;
calc_link = true;
plot_results = true;
if length(Cfg.Lat_vec)*length(Cfg.Lon_vec) > 8 % automate le decision
    use_parallel = true;
else
    use_parallel = false;
end

metrics = coverage_simulator_function(Cfg, plot_results,use_parallel,calc_link);
show_interactive = true;
save_fig = false;

show_constellation(Cfg, show_interactive, save_fig)
fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);
