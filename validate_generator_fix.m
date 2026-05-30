% VALIDATE_GENERATOR_FIX
% Validates the RAAN-budget seam-ratio fix on the three border altitudes
% where calculate_walker_star_inclined predicts a smaller constellation than
% SMAD, and where the old (ECA-ratio) generator produced a microscopic gap
% that made the gridsearch reject the inclined-formula constellation.
%
% Border altitudes (from the completed gridsearch sweep):
%   600 km  – SMAD T=132, inclined T=126, numerical found T=132 (gap ≈+0.007°)
%   760 km  – SMAD T=95,  inclined T=90,  numerical found T=95
%  1070 km  – SMAD T=60,  inclined T=56,  numerical found T=60
%
% For each altitude the script:
%   1. Computes [P,S,T] from calculate_walker_star_inclined(h, lat, elev, 90)
%      and calculate_walker_star(h, lat, elev) for comparison.
%   2. Runs the FULL detailed simulation — identical parameters to the
%      gridsearch detailed tier:  2150 UEs  |  350 h  |  660 s sampling
%      Lon range [-180,180] (global, matching gridsearch).
%   3. Reports worst_coverage_percent and PASS/FAIL vs 99.999%.
%
% PASS at every altitude means the fixed generator is consistent with the
% analytical formula, and the gridsearch will find T ≤ T_inclined.
%
% How to run:
%   matlab -nosplash -nodesktop -batch "validate_generator_fix" > validate_log.txt
%
% Branch: fix/raan-seam-ratio

clear; close all; clc;
addpath('functions');

%% ---- Configuration — must match Find_valid_walker_star_gridsearch.m -----
% Border altitudes only: these are the three cases that FAILED with the old
% ECA-ratio generator and should PASS with the RAAN-budget fix.
heights_km  = [600, 760, 1070];

lat         = 54 + 35/60;          % southern coverage boundary (deg N)
elev        = 20;                   % minimum service elevation (deg)
Lat_range   = [lat, 83 + 40/60];   % same as gridsearch
Lon_range   = [-180, 180];          % global — identical to gridsearch UE generation

% Detailed-tier parameters (gridsearch uses these for final confirmation)
NumUEs      = 2150;
duration_h  = 350;
sample_time = 660;                  % s per time step
threshold   = 99.999;              % % coverage required to declare VALID

%% ---- Randomised start time (same approach as gridsearch workers) ---------
rng(42);   % fixed seed for reproducibility
StartTime = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC') + hours((rand - 0.5) * 48);
StopTime  = StartTime + hours(duration_h);
fprintf('Simulation epoch: %s  →  %s\n', char(StartTime), char(StopTime));

%% ---- Baseline Cfg fields ------------------------------------------------
BaseCfg.Min_elevation_UE           = elev;
BaseCfg.Lat_range_deg              = Lat_range;
BaseCfg.SampleTime                 = sample_time;
BaseCfg.WalkerStar                 = true;
BaseCfg.StartTime                  = StartTime;
BaseCfg.StopTime                   = StopTime;
% Link-budget fields not needed for coverage-only check (calc_link=false)
BaseCfg.FRF                        = 1;
BaseCfg.RU                         = 1;
BaseCfg.Use_P618                   = false;
BaseCfg.Modified_shannon           = false;
BaseCfg.Simple_Atmospheric_Loss_dB = 0;
BaseCfg.Share_bandwidth            = false;
BaseCfg.Target_PFD_MHz             = -123;
BaseCfg.DL.Direction               = "DL";
BaseCfg.DL.f                       = 12e9;
BaseCfg.DL.B                       = 250e6;
BaseCfg.DL.NF                      = 5;
BaseCfg.DL.G_rx                    = 33;
BaseCfg.DL.Tx_type                 = "array";
BaseCfg.DL.Rx_type                 = "array";
BaseCfg.DL.G_tx                    = 30;
BaseCfg.DL.Max_P_tx_dBm            = 40;
BaseCfg.DL.Max_EIRP_dBm            = 70;
BaseCfg.DL.Max_EIRP_dBm_Hz         = 70 - 10*log10(250e6);
BaseCfg.DL.BeamGrid.num_beams      = 1;
BaseCfg.DL.BeamGrid.u_center       = 0;
BaseCfg.DL.BeamGrid.v_center       = 0;

%% ---- Validation loop ----------------------------------------------------
n = numel(heights_km);
results = struct();

sep = repmat('=', 1, 95);
fprintf('\n%s\n', sep);
fprintf('  GENERATOR FIX VALIDATION  |  lat=%.4f°  elev=%d°  UEs=%d  duration=%d h\n', ...
    lat, elev, NumUEs, duration_h);
fprintf('  Altitudes: %s km  (border cases from gridsearch)\n', ...
    strjoin(string(heights_km), ' / '));
fprintf('%s\n\n', sep);

hdr = sprintf('%-7s | %-19s | %-19s | %-11s | %-9s | %s', ...
    'Alt(km)', 'SMAD (P / S / T)', 'Inclined90 (P/S/T)', 'Cov%', 'T_saving', 'Result');
fprintf('%s\n%s\n', hdr, repmat('-', 1, numel(hdr)));

all_pass = true;
for k = 1:n
    h = heights_km(k);

    % Analytical formulas
    [p_smad, s_smad, t_smad] = calculate_walker_star(h, lat, elev);
    [p_incl, s_incl, t_incl] = calculate_walker_star_inclined(h, lat, elev, 90);

    % Fresh UE set per altitude (matches gridsearch worker behaviour)
    [UE_lats, UE_lons] = generate_equal_ish_area_UEs(Lat_range, Lon_range, NumUEs);

    % Build config for the inclined-formula (smaller) constellation
    Cfg                  = BaseCfg;
    Cfg.Orbit_height     = h * 1e3;      % m
    Cfg.Num_planes       = p_incl;
    Cfg.Sats_per_plane   = s_incl;
    Cfg.Total_sats       = t_incl;
    Cfg.Inclination      = 90;
    Cfg.Phasing          = p_incl / 2;
    Cfg.Flat_UE_array.Lats = UE_lats;
    Cfg.Flat_UE_array.Lons = UE_lons;

    fprintf('\n  [%d/%d] h=%d km | P=%d S=%d T=%d | %d UEs | %.0f h ...\n', ...
        k, n, h, p_incl, s_incl, t_incl, NumUEs, duration_h);

    t_start = tic;
    try
        % use_toolbox=false → fast_walker_star_ecef (fixed RAAN-budget seam ratio)
        % calc_link=false → coverage-only, no link budget
        % reset_cache=true → always use fresh computation
        metrics = constellation_simulator(Cfg, false, false, [], false, true);
        cov = metrics.worst_coverage_percent;
    catch ME
        cov = NaN;
        fprintf('  ERROR: %s\n', ME.message);
    end
    elapsed = toc(t_start);

    valid = ~isnan(cov) && cov >= threshold;
    if ~valid, all_pass = false; end

    results(k).h       = h;
    results(k).p_smad  = p_smad; results(k).s_smad = s_smad; results(k).t_smad = t_smad;
    results(k).p_incl  = p_incl; results(k).s_incl = s_incl; results(k).t_incl = t_incl;
    results(k).cov     = cov;
    results(k).valid   = valid;
    results(k).elapsed_s = elapsed;

    t_saving = t_smad - t_incl;   % satellites saved vs SMAD
    if valid; res = 'PASS'; else; res = 'FAIL'; end
    fprintf('%-7d | P=%2d S=%2d T=%3d      | P=%2d S=%2d T=%3d      | %-11.5f | %-9d | %s  (%.0fs)\n', ...
        h, p_smad, s_smad, t_smad, p_incl, s_incl, t_incl, cov, t_saving, res, elapsed);
end

fprintf('%s\n', repmat('-', 1, numel(hdr)));

%% ---- Summary -----------------------------------------------------------
n_pass = sum([results.valid]);
n_fail = n - n_pass;
fprintf('\n%s\n', sep);
if all_pass
    fprintf('  RESULT: ALL %d border altitudes PASSED (threshold %.4f%%)\n', n, threshold);
    fprintf('\n  Interpretation:\n');
    fprintf('  The fixed RAAN-budget seam ratio produces gap-free constellations\n');
    fprintf('  for every calculate_walker_star_inclined(h,lat,elev,90) prediction.\n');
    fprintf('  Running the gridsearch on this branch will find T <= T_inclined at\n');
    fprintf('  600 km, 760 km, and 1070 km — consistent with the analytical formula.\n');
else
    fprintf('  RESULT: %d/%d PASSED  |  %d FAILED\n', n_pass, n, n_fail);
    fprintf('\n  Failed altitudes:\n');
    for k = 1:n
        if ~results(k).valid
            fprintf('    %d km  P=%d S=%d T=%d  cov=%.5f%%\n', ...
                results(k).h, results(k).p_incl, results(k).s_incl, ...
                results(k).t_incl, results(k).cov);
        end
    end
end
fprintf('%s\n', sep);

%% ---- Save results -------------------------------------------------------
out_dir = fullfile('simulation_output', ...
    sprintf('validate_generator_fix_%s', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
mkdir(out_dir);
save(fullfile(out_dir, 'validation_results.mat'), 'results', 'heights_km', 'threshold', ...
    'NumUEs', 'duration_h', 'sample_time');
fprintf('\nResults saved to: %s\n', out_dir);
