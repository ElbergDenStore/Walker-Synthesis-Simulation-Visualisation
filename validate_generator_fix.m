% VALIDATE_GENERATOR_FIX
% Validates that the fixed asymmetrical_walker_star_generation (RAAN-budget
% seam ratio) produces constellations that the simulator accepts for EVERY
% altitude predicted by calculate_walker_star_inclined.
%
% For each altitude the script:
%   1. Computes the analytical minimum from calculate_walker_star_inclined(h,lat,elev,90)
%      and calculate_walker_star(h,lat,elev) for reference.
%   2. Runs a fast coverage simulation (use_SGP=false, 3-hour, 200 UEs) for
%      the inclined-analytical constellation.
%   3. Reports worst_coverage_percent and whether it meets the 99.999% threshold.
%
% PASS = analytical constellation is valid → numerical gridsearch can always
%        find a solution with T ≤ T_analytical (it just needs to test it).
%
% How to run:
%   matlab -nosplash -nodesktop -batch "validate_generator_fix" > validate_log.txt
%
% Branch: fix/raan-seam-ratio

clear; close all; clc;
addpath('functions');

%% ---- Configuration -------------------------------------------------------
heights_km   = 500:100:1200;   % altitudes to test
lat          = 54 + 35/60;     % Denmark minimum
elev         = 20;
Lat_range    = [lat, 83 + 40/60];
Lon_range    = [-73 - 10/60, 33 + 30/60];
NumUEs       = 200;
duration_h   = 3;              % short run — only checking for holes
sample_time  = 660;            % s per step
threshold    = 99.999;         % % coverage required to declare VALID

StartTime = datetime('2025-06-01 12:00:00', 'TimeZone', 'UTC');
StopTime  = StartTime + hours(duration_h);

%% ---- Generate UEs once (same set for all altitudes) ----------------------
[UE_lats, UE_lons] = generate_equal_ish_area_UEs(Lat_range, Lon_range, NumUEs);

%% ---- Baseline Cfg fields (no link budget — coverage only) ---------------
BaseCfg.Min_elevation_UE         = elev;
BaseCfg.Lat_range_deg            = Lat_range;
BaseCfg.SampleTime               = sample_time;
BaseCfg.StartTime                = StartTime;
BaseCfg.StopTime                 = StopTime;
BaseCfg.WalkerStar               = true;
BaseCfg.Flat_UE_array.Lats       = UE_lats;
BaseCfg.Flat_UE_array.Lons       = UE_lons;
% Link budget fields not needed for coverage-only check
BaseCfg.FRF                      = 1;
BaseCfg.RU                       = 1;
BaseCfg.Use_P618                 = false;
BaseCfg.Modified_shannon         = false;
BaseCfg.Simple_Atmospheric_Loss_dB = 0;
BaseCfg.Share_bandwidth          = false;
BaseCfg.Target_PFD_MHz           = -123;
BaseCfg.DL.Direction             = "DL";
BaseCfg.DL.f                     = 12e9;
BaseCfg.DL.B                     = 250e6;
BaseCfg.DL.NF                    = 5;
BaseCfg.DL.G_rx                  = 33;
BaseCfg.DL.Tx_type               = "array";
BaseCfg.DL.Rx_type               = "array";
BaseCfg.DL.G_tx                  = 30;
BaseCfg.DL.Max_P_tx_dBm          = 40;
BaseCfg.DL.Max_EIRP_dBm          = 70;
BaseCfg.DL.Max_EIRP_dBm_Hz       = 70 - 10*log10(250e6);
BaseCfg.DL.BeamGrid.num_beams    = 1;
BaseCfg.DL.BeamGrid.u_center     = 0;
BaseCfg.DL.BeamGrid.v_center     = 0;

%% ---- Run validation loop ------------------------------------------------
n = numel(heights_km);
results = struct();

sep = repmat('=', 1, 90);
fprintf('\n%s\n', sep);
fprintf('  GENERATOR FIX VALIDATION  |  lat=%.4f°  elev=%.1f°  UEs=%d  duration=%dh\n', ...
    lat, elev, NumUEs, duration_h);
fprintf('%s\n\n', sep);

hdr = sprintf('%-8s | %-18s | %-18s | %-10s | %-8s | %s', ...
    'Alt(km)', 'SMAD (P/S/T)', 'Inclined90 (P/S/T)', 'Cov%', 'T_margin', 'Result');
fprintf('%s\n%s\n', hdr, repmat('-', 1, numel(hdr)));

all_pass = true;
for k = 1:n
    h = heights_km(k);

    % Analytical formulas
    [p_smad,  s_smad,  t_smad]  = calculate_walker_star(h, lat, elev);
    [p_incl,  s_incl,  t_incl]  = calculate_walker_star_inclined(h, lat, elev, 90);

    % Build config for the inclined-formula constellation
    Cfg              = BaseCfg;
    Cfg.Orbit_height = h * 1e3;
    Cfg.Num_planes   = p_incl;
    Cfg.Sats_per_plane = s_incl;
    Cfg.Total_sats   = t_incl;
    Cfg.Inclination  = 90;
    Cfg.Phasing      = p_incl / 2;

    % Fast simulation (pure-math path, no toolbox)
    try
        metrics = fast_coverage_simulator_function(Cfg, true, false, false);
        cov = metrics.worst_coverage_percent;
    catch ME
        cov = NaN;
        fprintf('  ERROR at h=%d km: %s\n', h, ME.message);
    end

    valid = ~isnan(cov) && cov >= threshold;
    if ~valid, all_pass = false; end

    results(k).h       = h;
    results(k).p_smad  = p_smad;  results(k).s_smad = s_smad;  results(k).t_smad = t_smad;
    results(k).p_incl  = p_incl;  results(k).s_incl = s_incl;  results(k).t_incl = t_incl;
    results(k).cov     = cov;
    results(k).valid   = valid;

    t_margin = t_smad - t_incl;
    fprintf('%-8d | P=%2d S=%2d T=%3d     | P=%2d S=%2d T=%3d     | %-10.4f | %-8d | %s\n', ...
        h, p_smad, s_smad, t_smad, p_incl, s_incl, t_incl, cov, t_margin, ...
        char(valid * "PASS" + ~valid * "FAIL"));
end

fprintf('%s\n', repmat('-', 1, numel(hdr)));

%% ---- Summary -----------------------------------------------------------
n_pass = sum([results.valid]);
n_fail = n - n_pass;
fprintf('\n%s\n', sep);
if all_pass
    fprintf('  RESULT: ALL %d altitudes PASSED (%.4f%% threshold)\n', n, threshold);
    fprintf('  The fixed generator produces valid constellations for every\n');
    fprintf('  calculate_walker_star_inclined prediction at i=90°.\n');
    fprintf('  The numerical gridsearch is guaranteed to find T <= T_inclined.\n');
else
    fprintf('  RESULT: %d/%d PASSED  |  %d FAILED\n', n_pass, n, n_fail);
    fprintf('  Failed altitudes:\n');
    for k = 1:n
        if ~results(k).valid
            fprintf('    %d km  (P=%d S=%d T=%d  cov=%.4f%%)\n', ...
                results(k).h, results(k).p_incl, results(k).s_incl, ...
                results(k).t_incl, results(k).cov);
        end
    end
    fprintf('\n  Possible causes:\n');
    fprintf('  - Simulation duration too short (increase duration_h)\n');
    fprintf('  - Remaining issue in generator or formula\n');
end
fprintf('%s\n', sep);

%% ---- Save results -------------------------------------------------------
out_dir = fullfile('simulation_output', ...
    sprintf('validate_generator_fix_%s', char(datetime('now','Format','yyyyMMdd_HHmmss'))));
mkdir(out_dir);
save(fullfile(out_dir, 'validation_results.mat'), 'results', 'heights_km', 'threshold');
fprintf('\nResults saved to: %s\n', out_dir);
