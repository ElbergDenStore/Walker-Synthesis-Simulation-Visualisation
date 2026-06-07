% DIAGNOSE_COVERAGE_GAP
% Identifies the exact UE, time-step and orbital geometry behind the single
% coverage gap found by validate_generator_fix at 760 km and 1070 km.
%
% Reproduces the identical RNG state and UE positions, runs the full 350-h
% simulation in UE batches (avoids 8+ GB tensor), then calls show_constellation
% at the failure instant so you can inspect the geometry visually.
%
% Run from MATLAB command window (not -batch) so show_constellation can open
% an interactive satellite scenario viewer.

clear; close all; clc;
% Add the project to the MATLAB path (robust to the script's folder depth).
repo_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
addpath(fullfile(repo_root, 'functions'));
path_setup();

%% ---- Exact config from validate_generator_fix ---------------------------
rng(42);
StartTime  = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC') + hours((rand - 0.5) * 48);
lat        = 54 + 35/60;
elev       = 20;
Lat_range  = [lat, 83 + 40/60];
Lon_range  = [-180, 180];
NumUEs     = 2150;
duration_h = 350;
sample_time = 660;
StopTime   = StartTime + hours(duration_h);

fprintf('Epoch: %s\n\n', char(StartTime));

% UEs are deterministic (generate_equal_area_ues uses no random numbers)
[UE_lats, UE_lons] = generate_equal_area_ues(Lat_range, Lon_range, NumUEs);
ue_pos_ecef = lla2ecef([UE_lats, UE_lons, zeros(NumUEs, 1)]);

% Pre-build ENU rotation matrices [3 x 3 x NumUEs]
slat = sind(UE_lats'); clat = cosd(UE_lats');
slon = sind(UE_lons'); clon = cosd(UE_lons');
R_enu = zeros(3, 3, NumUEs);
R_enu(1,1,:) = -slon;       R_enu(1,2,:) = clon;        R_enu(1,3,:) = 0;
R_enu(2,1,:) = -slat.*clon; R_enu(2,2,:) = -slat.*slon; R_enu(2,3,:) = clat;
R_enu(3,1,:) = clat.*clon;  R_enu(3,2,:) = clat.*slon;  R_enu(3,3,:) = slat;

%% ---- Process each failing altitude ------------------------------------
heights_km = [760, 1070];
BATCH = 50;   % UEs processed at once — keeps memory under ~0.5 GB per batch

for h = heights_km
    [p_incl, s_incl, t_incl] = calculate_walker_star(h, lat, elev, 90);
    fprintf('=== %d km  P=%d S=%d T=%d ===\n', h, p_incl, s_incl, t_incl);

    % Build satellite positions using MATLAB two-body Keplerian — identical to
    % show_constellation, so visualiser and coverage computation use the same orbits.
    sc_diag = satelliteScenario;
    sc_diag.StartTime  = StartTime;
    sc_diag.StopTime   = StopTime;
    sc_diag.SampleTime = sample_time;
    sats_diag = generate_walker_star_scenario(sc_diag, h*1e3, 90, p_incl, s_incl, ...
                                                    elev, "two-body-keplerian", lat);
    [sat_pos_raw, ~, simTimes_diag] = states(sats_diag, "CoordinateFrame", "ECEF");
    % states() returns [3 x nT x numSats]; permute to [3 x numSats x nT]
    sat_pos  = permute(sat_pos_raw, [1, 3, 2]);
    time_steps = seconds(simTimes_diag - StartTime);
    nT         = length(time_steps);
    num_sats   = p_incl * s_incl;

    % Accumulate per-UE coverage count [NumUEs x nT]
    cov_counts = zeros(NumUEs, nT, 'uint16');

    n_batches = ceil(NumUEs / BATCH);
    for b = 1:n_batches
        idx = (b-1)*BATCH + 1 : min(b*BATCH, NumUEs);
        nb  = numel(idx);

        % [3 x num_sats x nT x nb]
        sat_4d = reshape(sat_pos,          3, num_sats, nT, 1);
        ue_4d  = reshape(ue_pos_ecef(idx,:)', 3, 1,       1,  nb);
        vec_4d = sat_4d - ue_4d;

        % Rotate to ENU for this batch
        vec_pages = reshape(vec_4d, 3, num_sats * nT, nb);
        enu_pages = pagemtimes(R_enu(:,:,idx), vec_pages);
        enu_4d    = reshape(enu_pages, 3, num_sats, nT, nb);

        U_4d  = squeeze(enu_4d(3, :, :, :));     % [num_sats x nT x nb]
        r_4d  = squeeze(sqrt(sum(enu_4d.^2, 1))); % [num_sats x nT x nb]
        el_4d = asind(U_4d ./ r_4d);              % [num_sats x nT x nb]

        cov_counts(idx, :) = uint16(squeeze(sum(el_4d >= elev, 1))'); % [nb x nT]
    end

    % Coverage fraction per UE
    cov_frac = sum(cov_counts > 0, 2) / nT * 100;  % [NumUEs x 1]
    [worst_cov, worst_ue] = min(cov_frac);

    % Time steps where that UE has no coverage
    fail_steps = find(cov_counts(worst_ue, :) == 0);

    fprintf('  Worst UE #%d\n', worst_ue);
    fprintf('  Location   : lat=%.5f°  lon=%.5f°\n', UE_lats(worst_ue), UE_lons(worst_ue));
    fprintf('  Coverage   : %.5f%%  (%d/%d steps covered)\n', worst_cov, nT - numel(fail_steps), nT);
    fprintf('  Gap steps  : %d\n', numel(fail_steps));

    for fi = 1:numel(fail_steps)
        fs         = fail_steps(fi);
        fail_time  = StartTime + seconds(time_steps(fs));

        % Max elevation of any satellite at this step for this UE
        ue_ecef  = ue_pos_ecef(worst_ue, :)';
        vec_sats = squeeze(sat_pos(:, :, fs)) - ue_ecef;  % [3 x num_sats]
        R_this   = squeeze(R_enu(:, :, worst_ue));
        enu_sats = R_this * vec_sats;                      % [3 x num_sats]
        el_sats  = asind(enu_sats(3,:) ./ sqrt(sum(enu_sats.^2, 1)));

        [max_el, best_sat] = max(el_sats);

        fprintf('\n  Gap #%d @ step %d  (%s)\n', fi, fs, char(fail_time));
        fprintf('    Nearest sat: #%d   elevation = %.5f°  (threshold %.1f°)\n', ...
                best_sat, max_el, elev);
        fprintf('    Deficit    : %.5f°  below threshold\n', elev - max_el);

        % ---- Show constellation at failure instant --------------------
        fprintf('\n  Calling show_constellation at failure instant...\n');
        window = seconds(sample_time);
        CfgViz.WalkerStar        = true;
        CfgViz.Orbit_height      = h * 1e3;
        CfgViz.Inclination       = 90;
        CfgViz.Num_planes        = p_incl;
        CfgViz.Sats_per_plane    = s_incl;
        CfgViz.Total_sats        = t_incl;
        CfgViz.Phasing           = p_incl / 2;
        CfgViz.Min_elevation_UE  = elev;
        CfgViz.Lat_range_deg     = Lat_range;
        CfgViz.SampleTime        = sample_time/10;
        % EpochTime must match the diagnostic scenario epoch so satellite positions
        % at fail_time are identical in both computations.
        CfgViz.EpochTime         = StartTime;
        CfgViz.StartTime         = fail_time;          % viewer jumps here on open
        CfgViz.StopTime          = fail_time + window;
        % Only the failing UE — beams are suppressed when NumUEs is large
        CfgViz.Flat_UE_array.Lats = UE_lats(worst_ue);
        CfgViz.Flat_UE_array.Lons = UE_lons(worst_ue);
        show_constellation(CfgViz, true, false);
    end
    fprintf('\n');
end
