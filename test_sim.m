try
    addpath('functions');
    Cfg = get_cfg(1000,'walkerdelta','small','medium');
    calc_link = true; 
    plot_results = true; 
    use_parallel = false; 
    metrics = coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link);
    disp('SIMULATION SUCCESS!');
catch ME
    disp(getReport(ME));
end
quit;
