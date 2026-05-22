% compare_toolbox_vs_math.m
% Compares MATLAB two-body-keplerian (toolbox) against analytical two-body
% math. Works for BOTH WalkerStar and Walker Delta — the orbital mechanics
% (circular two-body, ECI→ECEF via GMST) are identical regardless.
clc; clear; close all;
addpath('functions');

%% ── Configuration ────────────────────────────────────────────────────────
Cfg.WalkerStar       = true;   % true = Walker Star,  false = Walker Delta
Cfg.Orbit_height     = 550e3;  % metres
Cfg.Inclination      = 90;
Cfg.Num_planes       = 7;
Cfg.Sats_per_plane   = 23;
Cfg.Total_sats       = Cfg.Num_planes * Cfg.Sats_per_plane;
Cfg.Phasing          = Cfg.Num_planes / 2;  % only used for Walker Delta
Cfg.Min_elevation_UE = 20;
Cfg.SampleTime       = 60;
Cfg.Lat_range_deg    = [54+(35/60), 83+(40/60)];
StartTime            = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
Cfg.StartTime        = StartTime;

r_earth = 6378.14e3;          % metres – matches toolbox paths throughout
a_m     = r_earth + Cfg.Orbit_height;
mu      = 3.986004418e14;
inc     = deg2rad(Cfg.Inclination);
P = Cfg.Num_planes;  S = Cfg.Sats_per_plane;  T = P * S;

%% ── Build orbital elements (RAAN, ν₀) for every satellite ───────────────
% The two-body math is the same for Walker Star and Walker Delta; only the
% initial element assignment differs.
raans = zeros(1, T);
nu0s  = zeros(1, T);

if Cfg.WalkerStar
    % Replicate asymmetrical_walker_star_generation element assignment
    orbit_height_km = Cfg.Orbit_height / 1000;
    Re_eq_km = 6378.14;
    Rs_km    = Re_eq_km + orbit_height_km;
    a_wgs = 6378.137; b_wgs = 6356.7523142;
    lat_r = deg2rad(Cfg.Lat_range_deg(1));
    Re_km = sqrt((a_wgs^4*cos(lat_r)^2 + b_wgs^4*sin(lat_r)^2) / ...
                 (a_wgs^2*cos(lat_r)^2 + b_wgs^2*sin(lat_r)^2));
    alpha        = asind((Re_km / Rs_km) * cosd(Cfg.Min_elevation_UE));
    lambda_max   = deg2rad(180 - (90 + Cfg.Min_elevation_UE + alpha));
    S_rad        = 2*pi / S;
    lambda_street = acos(min(1, cos(lambda_max) / cos(S_rad / 2)));
    seam_ratio   = (2*lambda_street) / (lambda_street + lambda_max);
    co_rot_spacing = 180 / (P - 1 + seam_ratio);   % degrees
    in_plane_spc   = 360 / S;                        % degrees
    phase_shift    = in_plane_spc / 2;               % degrees
    k = 1;
    for p = 1:P
        for s = 1:S
            raans(k) = deg2rad((p-1) * co_rot_spacing);
            nu0s(k)  = deg2rad(mod((s-1)*in_plane_spc + (p-1)*phase_shift, 360));
            k = k + 1;
        end
    end
else
    % Walker Delta element assignment  (matches walkerDelta() and fast_walker_ecef)
    F = Cfg.Phasing;
    k = 1;
    for p = 1:P
        for s = 1:S
            raans(k) = (p-1) * 2*pi/P;
            nu0s(k)  = deg2rad(((s-1)*360/S + (p-1)*F*360/T));
            k = k + 1;
        end
    end
end

%% ── TEST 1: Direct satellite position comparison at t = 0 ────────────────
fprintf('\n--- TEST 1: Satellite position error at t=0 ---\n');

% GMST at StartTime
JD       = juliandate(StartTime);
theta_g0 = deg2rad(mod(280.46061837 + 360.98564736629*(JD - 2451545.0), 360));

% Analytical (two-body math) positions
pos_math = zeros(3, T);
for k = 1:T
    RAAN  = raans(k);   nu = nu0s(k);
    x_orb = a_m * cos(nu);   y_orb = a_m * sin(nu);
    X_eci =  x_orb*cos(RAAN) - y_orb*cos(inc)*sin(RAAN);
    Y_eci =  x_orb*sin(RAAN) + y_orb*cos(inc)*cos(RAAN);
    Z_eci =  y_orb * sin(inc);
    pos_math(1,k) =  X_eci*cos(theta_g0) + Y_eci*sin(theta_g0);
    pos_math(2,k) = -X_eci*sin(theta_g0) + Y_eci*cos(theta_g0);
    pos_math(3,k) =  Z_eci;
end

% MATLAB toolbox positions
sc_test            = satelliteScenario;
sc_test.StartTime  = StartTime;
sc_test.StopTime   = StartTime + seconds(Cfg.SampleTime);
sc_test.SampleTime = Cfg.SampleTime;
if Cfg.WalkerStar
    sats_tb = asymmetrical_walker_star_generation(sc_test, Cfg.Orbit_height, ...
        Cfg.Inclination, P, S, Cfg.Min_elevation_UE, "two-body-keplerian", Cfg.Lat_range_deg(1));
else
    sats_tb = walkerDelta(sc_test, a_m, Cfg.Inclination, T, P, Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="two-body-keplerian");
end
pos_raw     = states(sats_tb, "CoordinateFrame", "ECEF");  % [3 x nT x nSats]
pos_toolbox = squeeze(pos_raw(:, 1, :));                   % [3 x nSats] at t=0

% Nearest-neighbour error (handles any satellite ordering difference)
pos_err = zeros(1, T);
for k = 1:T
    pos_err(k) = min(vecnorm(pos_toolbox - pos_math(:,k), 2, 1));
end
fprintf('Max position error at t=0:  %.2f m\n', max(pos_err));
fprintf('Mean position error at t=0: %.2f m\n', mean(pos_err));

%% ── TEST 2 & 3: Coverage statistics ──────────────────────────────────────
NumUEs = 200;
[UE_lats, UE_lons] = generate_equal_ish_area_UEs(Cfg.Lat_range_deg, [-180,180], NumUEs);
Cfg.Flat_UE_array.Lats = UE_lats;
Cfg.Flat_UE_array.Lons = UE_lons;

if Cfg.WalkerStar
    fprintf('\n[NOTE] WalkerStar=true: fast_coverage_simulator_function falls back\n');
    fprintf('       to the toolbox for use_SGP=false. Coverage diff will be ~0.\n');
end

fprintf('\n--- TEST 2: Coverage statistics (3 Hours) ---\n');
Cfg.StopTime = StartTime + hours(3);
m2_toolbox = fast_coverage_simulator_function(Cfg, true,  false, true);
m2_math    = fast_coverage_simulator_function(Cfg, false, false, false);
fprintf('Worst Coverage:\n  Toolbox: %.4f%%\n  Math:    %.4f%%\n', ...
    m2_toolbox.worst_coverage_percent, m2_math.worst_coverage_percent);

fprintf('\n--- TEST 3: Coverage statistics (24 Hours) ---\n');
Cfg.StopTime = StartTime + hours(24);
m3_toolbox = fast_coverage_simulator_function(Cfg, true,  false, true);
m3_math    = fast_coverage_simulator_function(Cfg, false, false, false);
fprintf('Worst Coverage:\n  Toolbox: %.4f%%\n  Math:    %.4f%%\n', ...
    m3_toolbox.worst_coverage_percent, m3_math.worst_coverage_percent);

%% ── Summary ──────────────────────────────────────────────────────────────
fprintf('\n--- Summary ---\n');
if     max(pos_err) <   100, fprintf('[POS]     VERY close at t=0 (< 100 m).\n');
elseif max(pos_err) <  1000, fprintf('[POS]     Close at t=0 (< 1 km).\n');
else,  fprintf('[POS]     Noticeable offset at t=0 (%.2f km).\n', max(pos_err)/1e3);
end
fprintf('[3H  COV] Difference: %.4f%%\n', abs(m2_toolbox.worst_coverage_percent - m2_math.worst_coverage_percent));
fprintf('[24H COV] Difference: %.4f%%\n', abs(m3_toolbox.worst_coverage_percent - m3_math.worst_coverage_percent));
