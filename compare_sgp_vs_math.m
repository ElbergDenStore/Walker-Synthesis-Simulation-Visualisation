% compare_sgp_vs_math.m
clc; clear; close all;

% Configuration
Cfg.Orbit_height = 550e3;
Cfg.Inclination = 90;
Cfg.Num_planes = 7;
Cfg.Sats_per_plane = 23;
Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
Cfg.Phasing = Cfg.Num_planes/2;
Cfg.WalkerStar = true;
Cfg.Min_elevation_UE = 20;
Cfg.SampleTime = 60;
Lat_range_deg = [54+(35/60), 83+(40/60)];
% Lon_range_deg = [-(73+(10/60)), 33+(30/60)];
Lon_range_deg = [-180, 180];
StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
StopTime  = datetime('3-Jun-2025 12:59:59', 'TimeZone', 'UTC'); % 48 hours



NumUEs = 200;

[UE_lats, UE_lons] = generate_equal_ish_area_UEs(Lat_range_deg, Lon_range_deg, NumUEs);
Cfg.Flat_UE_array.Lats       = UE_lats;
Cfg.Flat_UE_array.Lons       = UE_lons;
Cfg.StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');

%% TEST 1: Short-term accuracy (10 Minutes)
fprintf('\n--- TEST 1: Short-term accuracy (10 Minutes) ---\n');
Cfg.StopTime = Cfg.StartTime + minutes(10);
m1_sgp  = fast_coverage_simulator_function(Cfg, true, false, true);
m1_math = fast_coverage_simulator_function(Cfg, true, false, false);

% Max Range Diff
max_r = 0;
for i = 1:length(m1_sgp.SimData)
    v = m1_sgp.SimData(i).Num_visible > 0 & m1_math.SimData(i).Num_visible > 0;
    if any(v)
        max_r = max(max_r, max(abs(m1_sgp.SimData(i).Range(v) - m1_math.SimData(i).Range(v))));
    end
end
fprintf(sprintf("start difference %.2f meters \n",abs(m1_sgp.SimData(1).Range(1) - m1_math.SimData(1).Range(1))))
fprintf('Initial Max Range Difference: %.2f meters\n', max_r);


%% TEST 2: Long-term statistical (24 Hours)
fprintf('\n--- TEST 2: Long-term statistical (24 Hours) ---\n');
Cfg.StopTime = Cfg.StartTime + hours(24);
m2_sgp  = fast_coverage_simulator_function(Cfg, true, false, false, true);
m2_math = fast_coverage_simulator_function(Cfg, true, false, false, false);

fprintf('Worst Coverage Percent:\n  SGP:  %.4f%%\n  MATH: %.4f%%\n', m2_sgp.worst_coverage_percent, m2_math.worst_coverage_percent);
fprintf('Worst Gap Minutes:\n  SGP:  %.2f min\n  MATH: %.2f min\n', m2_sgp.worst_gap_minutes, m2_math.worst_gap_minutes);

fprintf('\nSummary of findings:\n');
if max_r < 100
    fprintf('- [10 MIN] Initial math is VERY close to SGP (under 100m error).\n');
elseif max_r < 1000
    fprintf('- [10 MIN] Initial math is close (under 1km error).\n');
else
    fprintf('- [10 MIN] Initial math has noticeable drift even at start (%.2f km error).\n', max_r/1e3);
end

coverage_diff = abs(m2_sgp.worst_coverage_percent - m2_math.worst_coverage_percent);
if coverage_diff < 0.5
    fprintf('- [24 HOUR] Coverage statistics are very similar (within 0.5%%).\n');
elseif coverage_diff < 2
    fprintf('- [24 HOUR] Coverage statistics are fairly close (within 2%%).\n');
else
    fprintf('- [24 HOUR] Coverage statistics diverge significantly (%.2f%% difference).\n', coverage_diff);
end
