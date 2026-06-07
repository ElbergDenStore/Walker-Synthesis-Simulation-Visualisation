% validate_propagator_accuracy.m
% Validates the pure-math ECEF generator (fast_walker_ecef, Walker Star + Delta)
% against the MATLAB Satellite Toolbox two-body-keplerian propagator.
% Runs three tests for BOTH Walker Star and Walker Delta:
%   TEST 1 – Direct satellite position error at t = 0
%   TEST 2 – Worst-case coverage difference over 3 hours
%   TEST 3 – Worst-case coverage difference over 24 hours
clc; clear; close all;
% Add the project to the MATLAB path (robust to the script's folder depth).
repo_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
addpath(fullfile(repo_root, 'functions'));
path_setup();

StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
NumUEs    = 200;

%% ── Constellation Configurations ─────────────────────────────────────────
cfgStar.WalkerStar       = true;
cfgStar.Orbit_height     = 510e3;     % m
cfgStar.Inclination      = 90;
cfgStar.Num_planes       = 7;
cfgStar.Sats_per_plane   = 22;
cfgStar.Total_sats       = cfgStar.Num_planes * cfgStar.Sats_per_plane;
cfgStar.Phasing          = 0;         % unused for Walker Star
cfgStar.Min_elevation_UE = 20;
cfgStar.SampleTime       = 60;
cfgStar.Lat_range_deg    = [54+(35/60), 83+(40/60)];
cfgStar.StartTime        = StartTime;

cfgDelta.WalkerStar       = false;
cfgDelta.Orbit_height     = 1000e3;    % m
cfgDelta.Inclination      = 75;
cfgDelta.Num_planes       = 4;
cfgDelta.Sats_per_plane   = 13;
cfgDelta.Total_sats       = cfgDelta.Num_planes * cfgDelta.Sats_per_plane;
cfgDelta.Phasing          = 2;
cfgDelta.Min_elevation_UE = 20;
cfgDelta.SampleTime       = 60;
cfgDelta.Lat_range_deg    = [54+(35/60), 83+(40/60)];
cfgDelta.StartTime        = StartTime;

configs = {cfgStar, cfgDelta};
labels  = {'Walker Star', 'Walker Delta'};
summary = cell(numel(configs), 1);

%% ── Run all three tests for each constellation type ─────────────────────
for ci = 1:numel(configs)
    Cfg = configs{ci};
    fprintf('\n%s\n', repmat('=', 1, 65));
    fprintf('  %s  (%d sats, h=%g km, i=%g deg)\n', ...
        labels{ci}, Cfg.Total_sats, Cfg.Orbit_height/1e3, Cfg.Inclination);
    fprintf('%s\n', repmat('=', 1, 65));

    r_earth = 6378.137e3;
    a_m  = r_earth + Cfg.Orbit_height;
    inc  = deg2rad(Cfg.Inclination);
    P = Cfg.Num_planes;  S = Cfg.Sats_per_plane;  T = P * S;

    %% TEST 1: Direct satellite position comparison at t = 0 ---------------
    fprintf('\n--- TEST 1: Satellite position error at t=0 ---\n');

    [raans, nu0s] = build_elements(Cfg);

    JD       = juliandate(StartTime);
    theta_g0 = deg2rad(mod(280.46061837 + 360.98564736629*(JD - 2451545.0), 360));

    % Analytical (pure-math) positions
    pos_math = zeros(3, T);
    for k = 1:T
        RAAN  = raans(k);  nu = nu0s(k);
        x_orb = a_m * cos(nu);   y_orb = a_m * sin(nu);
        X_eci =  x_orb*cos(RAAN) - y_orb*cos(inc)*sin(RAAN);
        Y_eci =  x_orb*sin(RAAN) + y_orb*cos(inc)*cos(RAAN);
        Z_eci =  y_orb * sin(inc);
        pos_math(1,k) =  X_eci*cos(theta_g0) + Y_eci*sin(theta_g0);
        pos_math(2,k) = -X_eci*sin(theta_g0) + Y_eci*cos(theta_g0);
        pos_math(3,k) =  Z_eci;
    end

    % Toolbox positions
    sc_test            = satelliteScenario;
    sc_test.StartTime  = StartTime;
    sc_test.StopTime   = StartTime + seconds(Cfg.SampleTime);
    sc_test.SampleTime = Cfg.SampleTime;
    if Cfg.WalkerStar
        sats_tb = generate_walker_star_scenario(sc_test, Cfg.Orbit_height, ...
            Cfg.Inclination, P, S, Cfg.Min_elevation_UE, "two-body-keplerian", ...
            Cfg.Lat_range_deg(1));
    else
        sats_tb = walkerDelta(sc_test, a_m, Cfg.Inclination, T, P, Cfg.Phasing, ...
            Name="S4D", OrbitPropagator="two-body-keplerian");
    end
    pos_raw     = states(sats_tb, "CoordinateFrame", "ECEF");  % [3 x nT x nSats]
    pos_toolbox = squeeze(pos_raw(:, 1, :));                   % [3 x nSats] at t=0

    % Nearest-neighbour matching (tolerates any satellite ordering difference)
    pos_err = zeros(1, T);
    for k = 1:T
        pos_err(k) = min(vecnorm(pos_toolbox - pos_math(:,k), 2, 1));
    end
    fprintf('Max position error at t=0:  %.2f m\n', max(pos_err));
    fprintf('Mean position error at t=0: %.2f m\n', mean(pos_err));

    %% TEST 2 & 3: Coverage statistics -------------------------------------
    [UE_lats, UE_lons] = generate_equal_area_ues(Cfg.Lat_range_deg, [-180, 180], NumUEs);
    Cfg.Flat_UE_array.Lats = UE_lats;
    Cfg.Flat_UE_array.Lons = UE_lons;

    fprintf('\n--- TEST 2: Coverage statistics (3 Hours) ---\n');
    Cfg.StopTime = StartTime + hours(3);
    m2_toolbox = Constellation_simulator(Cfg, false, false, true);
    m2_math    = Constellation_simulator(Cfg, false, false, false);
    fprintf('Worst Coverage:\n  Toolbox: %.4f%%\n  Math:    %.4f%%\n', ...
        m2_toolbox.worst_coverage_percent, m2_math.worst_coverage_percent);

    fprintf('\n--- TEST 3: Coverage statistics (24 Hours) ---\n');
    Cfg.StopTime = StartTime + hours(24);
    m3_toolbox = Constellation_simulator(Cfg, false, false, true);
    m3_math    = Constellation_simulator(Cfg, false, false, false);
    fprintf('Worst Coverage:\n  Toolbox: %.4f%%\n  Math:    %.4f%%\n', ...
        m3_toolbox.worst_coverage_percent, m3_math.worst_coverage_percent);

    summary{ci} = struct( ...
        'label',       labels{ci}, ...
        'max_pos_err', max(pos_err), ...
        'diff_3h',     abs(m2_toolbox.worst_coverage_percent - m2_math.worst_coverage_percent), ...
        'diff_24h',    abs(m3_toolbox.worst_coverage_percent - m3_math.worst_coverage_percent));
end

%% ── Summary table ────────────────────────────────────────────────────────
fprintf('\n%s\n', repmat('=', 1, 65));
fprintf('  SUMMARY\n');
fprintf('%s\n', repmat('=', 1, 65));
fprintf('%-16s | %14s | %16s | %16s\n', ...
    'Constellation', 'Pos err (m)', '3h cov diff (%)', '24h cov diff (%)');
fprintf('%s\n', repmat('-', 1, 65));
for ci = 1:numel(summary)
    s = summary{ci};
    if     s.max_pos_err <   100, pos_tag = '  OK (<100 m)  ';
    elseif s.max_pos_err <  1000, pos_tag = '  OK (<1 km)   ';
    else,                         pos_tag = sprintf(' WARN (%.1f km)', s.max_pos_err/1e3);
    end
    fprintf('%-16s | %s | %16.4f | %16.4f\n', s.label, pos_tag, s.diff_3h, s.diff_24h);
end
fprintf('%s\n', repmat('=', 1, 65));
fprintf('Expectation: pos err <100 m, coverage diffs ~0%%.\n\n');

%% ── Helper: build initial orbital elements (RAAN + ν₀) ──────────────────
function [raans, nu0s] = build_elements(Cfg)
    P = Cfg.Num_planes;  S = Cfg.Sats_per_plane;  T = P * S;
    raans = zeros(1, T);
    nu0s  = zeros(1, T);

    if Cfg.WalkerStar
        % Replicate generate_walker_star_scenario element assignment
        orbit_height_km = Cfg.Orbit_height / 1000;
        Re_eq_km = 6378.137;
        Rs_km    = Re_eq_km + orbit_height_km;
        a_wgs = 6378.137;  b_wgs = 6356.7523142;
        lat_r = deg2rad(Cfg.Lat_range_deg(1));
        Re_km = sqrt((a_wgs^4*cos(lat_r)^2 + b_wgs^4*sin(lat_r)^2) / ...
                     (a_wgs^2*cos(lat_r)^2 + b_wgs^2*sin(lat_r)^2));
        alpha         = asind((Re_km / Rs_km) * cosd(Cfg.Min_elevation_UE));
        lambda_max    = deg2rad(180 - (90 + Cfg.Min_elevation_UE + alpha));
        S_rad         = 2*pi / S;
        lambda_street = acos(min(1, cos(lambda_max) / cos(S_rad / 2)));
        seam_ratio    = (2*lambda_street) / (lambda_street + lambda_max);
        co_rot_spacing = 180 / (P - 1 + seam_ratio);
        in_plane_spc   = 360 / S;
        phase_shift    = in_plane_spc / 2;
        k = 1;
        for p = 1:P
            for s = 1:S
                raans(k) = deg2rad((p-1) * co_rot_spacing);
                nu0s(k)  = deg2rad(mod((s-1)*in_plane_spc + (p-1)*phase_shift, 360));
                k = k + 1;
            end
        end
    else
        % Walker Delta element assignment (matches walkerDelta() and fast_walker_ecef)
        F = Cfg.Phasing;
        k = 1;
        for p = 1:P
            for s = 1:S
                raans(k) = (p-1) * 2*pi/P;
                nu0s(k)  = deg2rad((s-1)*360/S + (p-1)*F*360/T);
                k = k + 1;
            end
        end
    end
end
