% validate_simulators.m
clc; clear; close all;

fprintf('--- Setting up Validation Configuration ---\n');
Cfg.Orbit_height = 550e3;
Cfg.Inclination = 75;
Cfg.Num_planes = 8;
Cfg.Sats_per_plane = 8;
Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
Cfg.Phasing = 1;
Cfg.WalkerStar = false;
Cfg.Min_elevation_UE = 20;

Cfg.StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
Cfg.StopTime = datetime('2-Jun-2025 12:30:00', 'TimeZone', 'UTC'); 
Cfg.SampleTime = 60;

Cfg.Lat_vec = linspace(50, 60, 4);
Cfg.Lon_vec = linspace(5, 15, 4);
Cfg.Equal_UE_area = false;
Cfg.Accept_Flat_UE_array = false;

fprintf('Number of UEs roughly: %d\n', length(Cfg.Lat_vec)*length(Cfg.Lon_vec));

fprintf('\n--- Running OLD Simulator (coverage_simulator_function) ---\n');
tic;
old_metrics = coverage_simulator_function(Cfg, false, false, false);
old_time = toc;
fprintf('OLD Simulator took %.2f seconds.\n', old_time);

fprintf('\n--- Running NEW FAST Simulator (fast_coverage_simulator_function) ---\n');
tic;
new_metrics = fast_coverage_simulator_function(Cfg, true, false, false);
new_time = toc;
fprintf('NEW Simulator took %.2f seconds.\n', new_time);

fprintf('\n--- Comparing Results ---\n');
fprintf('Speedup: %.2fx faster!\n', old_time / new_time);

% Compare high level stats
fprintf('\nWorst Coverage Percent:\n  OLD: %.4f%%\n  NEW: %.4f%%\n', old_metrics.worst_coverage_percent, new_metrics.worst_coverage_percent);

disp(size(old_metrics.Num_visible));
disp(size(new_metrics.Num_visible));
diff_visible = max(abs(old_metrics.Num_visible(:) - new_metrics.Num_visible(:)));
fprintf('Max difference in Num_visible: %d satellites\n', diff_visible);

% Compare detailed SimData
max_range_diff = 0;
max_el_diff = 0;
max_az_diff = 0;

num_service_mismatches = 0;

for i = 1:length(old_metrics.SimData)
    % Compare ranges when has_service
    has_service_old = old_metrics.SimData(i).Num_visible > 0;
    has_service_new = new_metrics.SimData(i).Num_visible > 0;
    
    if any(has_service_old ~= has_service_new)
        num_service_mismatches = num_service_mismatches + sum(has_service_old ~= has_service_new);
    end
    
    valid_idx = has_service_old & has_service_new;
    
    if any(valid_idx)
        r_diff = max(abs(old_metrics.SimData(i).Range(valid_idx) - new_metrics.SimData(i).Range(valid_idx)));
        e_diff = max(abs(old_metrics.SimData(i).Elevation_deg(valid_idx) - new_metrics.SimData(i).Elevation_deg(valid_idx)));
        
        % Azimuth difference needs modulo 360 wrap logic (e.g., 359 vs 1 degree is a 2 degree difference)
        a_old = old_metrics.SimData(i).Azimuth_deg(valid_idx);
        a_new = new_metrics.SimData(i).Azimuth_deg(valid_idx);
        a_diff_raw = abs(a_old - a_new);
        a_diff = max(min(a_diff_raw, 360 - a_diff_raw));
        
        max_range_diff = max(max_range_diff, r_diff);
        max_el_diff = max(max_el_diff, e_diff);
        max_az_diff = max(max_az_diff, a_diff);
    end
end

fprintf('\nService Logic:\n');
fprintf('  Mismatching service timestamps: %d\n', num_service_mismatches);

fprintf('\nMax Absolute Numerical Differences (During Valid Service):\n');
fprintf('  Range:     %.6f meters\n', max_range_diff);
fprintf('  Elevation: %.6f degrees\n', max_el_diff);
fprintf('  Azimuth:   %.6f degrees\n', max_az_diff);

% Thresholds for "identical" 
% Range diffs < 100 meters due to SGP4 numerical differences vs ECEF
% Elev/Azim diffs < 0.1 deg
if num_service_mismatches == 0 && max_range_diff < 100 && max_el_diff < 0.05 && max_az_diff < 0.05
    fprintf('\nSUCCESS! The new simulator produces practically identical metrics to the original MATLAB AER.\n');
else
    fprintf('\nWARNING! Significant differences detected. MathWorks aer() uses WGS84 ellipsoid whereas fast pure-math may use perfect spheres, or there is a bug.\n');
end
