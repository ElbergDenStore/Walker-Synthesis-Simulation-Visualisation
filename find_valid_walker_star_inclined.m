clear; close all; clc;

% Ensure functions/ subdirectory is on the MATLAB path
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

%% Parameters
heights_km       = 500:0.1:1200;
Min_latitude_deg = 54 + 35/60;   % 54°35'N — southernmost Denmark
Min_elevation_UE = 20;            % degrees minimum elevation angle
inc_array        = 90:-1:84;      % 90°, 89°, 88°, 87°, 86°, 85°, 84°

n_heights = numel(heights_km);
n_inc     = numel(inc_array);

%% Compute results
total_sats_all = zeros(n_heights, n_inc);
fprintf('Computing Walker Star coverage (%d inclinations x %d altitudes)...\n', ...
    n_inc, n_heights);

for ii = 1:n_inc
    for hi = 1:n_heights
        [~, ~, ts] = calculate_walker_star_inclined( ...
            heights_km(hi), Min_latitude_deg, Min_elevation_UE, inc_array(ii));
        if isinf(ts); ts = NaN; end
        total_sats_all(hi, ii) = ts;
    end
    fprintf('  i = %d deg done\n', inc_array(ii));
end

%% Save results
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir  = fullfile(fileparts(mfilename('fullpath')), 'simulation_output', ...
                    ['Walker-Star-Inclined_' date_str]);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save(fullfile(out_dir, 'total_sats_all.mat'), 'total_sats_all', 'heights_km', 'inc_array', ...
     'Min_latitude_deg', 'Min_elevation_UE');
fprintf('Results saved to: %s\n', out_dir);
