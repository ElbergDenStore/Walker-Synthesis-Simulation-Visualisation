function plot_altitude_vs_latitude()
% PLOT_ALTITUDE_VS_LATITUDE
%   Shows how geodetic altitude (WGS84) varies with latitude for a circular
%   orbit at 1000 km nominal height.  Two inclinations are overlaid:
%     87° – Walker Star design inclination
%     90° – True polar
%
%   Panel A: full orbital period  (latitude traces –max → 0 → +max → 0 → –max)
%   Panel B: ascending quarter only (equator → peak latitude)
%
%   The variation is purely geometric (WGS84 oblateness).  No SGP4 / J2
%   perturbations are included so the signal is unambiguous.

HEIGHT_KM    = 1000;
r_earth      = 6378.137e3;          % WGS84 semi-major axis [m]
a            = r_earth + HEIGHT_KM*1e3;
inclinations = [90];
colors       = {'#0072BD'};
labels       = {'i = 90° (polar)'};
n_pts        = 50000;

u_full = linspace(0, 2*pi, n_pts);  % argument of latitude, full orbit

%% Pre-compute for each inclination
results = struct();
for k = 1:numel(inclinations)
    inc = deg2rad(inclinations(k));

    % ECI position (RAAN = 0; latitude/altitude independent of RAAN)
    x = a * cos(u_full);
    y = a * sin(u_full) .* cos(inc);
    z = a * sin(u_full) .* sin(inc);

    lla     = ecef2lla([x(:), y(:), z(:)]);   % WGS84 geodetic
    lat     = lla(:,1);                        % degrees
    alt_km  = lla(:,3) / 1e3;

    % Nominal (spherical) altitude reference
    r_geocentric = sqrt(x.^2 + y.^2 + z.^2);
    nominal_alt_km = (r_geocentric(:) - r_earth) / 1e3;

    results(k).lat            = lat;
    results(k).alt_km         = alt_km;
    results(k).nominal_alt_km = nominal_alt_km;
    results(k).deviation_km   = alt_km - nominal_alt_km;   % WGS84 effect

    % Ascending quarter: u in [0, pi/2]  → latitude 0° → peak
    asc_mask = u_full <= pi/2;
    results(k).lat_asc = lat(asc_mask);
    results(k).alt_asc = alt_km(asc_mask);
    results(k).dev_asc = results(k).deviation_km(asc_mask);

    [peak_lat, ~] = max(abs(lat));
    fprintf('i=%d°  |  peak geodetic lat: %.2f°  |  alt range: [%.4f, %.4f] km  |  WGS84 deviation: [%.4f, %.4f] km\n', ...
        inclinations(k), peak_lat, min(alt_km), max(alt_km), ...
        min(results(k).deviation_km), max(results(k).deviation_km));
end

set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontSize', 13);

%% ── Figure 1: Full orbital period ───────────────────────────────────────
f1 = figure('Name', 'Altitude vs Latitude – Full Orbit', ...
    'Color', 'w', 'Position', [60 100 1100 480]);

subplot(1,2,1);  % Geodetic altitude
hold on;
for k = 1:numel(inclinations)
    % Sort by latitude for a clean closed curve
    [lat_s, si] = sort(results(k).lat);
    alt_s = results(k).alt_km(si);
    plot(lat_s, alt_s, '-', 'Color', colors{k}, 'LineWidth', 1.8, 'DisplayName', labels{k});
end
yline(HEIGHT_KM, 'k--', 'LineWidth', 1.2, 'DisplayName', sprintf('Nominal %d km (spherical)', HEIGHT_KM));
xlabel('Geodetic Latitude (°)');
ylabel('Geodetic Altitude (km)');
title('Geodetic altitude – full orbit');
legend('Location', 'south');
grid on;

subplot(1,2,2);  % WGS84 deviation only
hold on;
for k = 1:numel(inclinations)
    [lat_s, si] = sort(results(k).lat);
    dev_s = results(k).deviation_km(si);
    plot(lat_s, dev_s, '-', 'Color', colors{k}, 'LineWidth', 1.8, 'DisplayName', labels{k});
end
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('Geodetic Latitude (°)');
ylabel('Geodetic alt − spherical alt (km)');
title('WGS84 deviation from spherical nominal');
legend('Location', 'south');
grid on;

%% ── Figure 2: Ascending pass only (equator → peak latitude) ─────────────
f2 = figure('Name', 'Altitude vs Latitude – Ascending Pass', ...
    'Color', 'w', 'Position', [60 580 1100 480]);

subplot(1,2,1);  % Geodetic altitude
hold on;
for k = 1:numel(inclinations)
    [lat_s, si] = sort(results(k).lat_asc);
    alt_s = results(k).alt_asc(si);
    plot(lat_s, alt_s, '-', 'Color', colors{k}, 'LineWidth', 1.8, 'DisplayName', labels{k});
end
yline(HEIGHT_KM, 'k--', 'LineWidth', 1.2, 'DisplayName', sprintf('Nominal %d km', HEIGHT_KM));
xlabel('Geodetic Latitude (°)');
ylabel('Geodetic Altitude (km)');
title('Geodetic altitude – ascending pass');
legend('Location', 'best');
grid on;

subplot(1,2,2);  % WGS84 deviation only
hold on;
for k = 1:numel(inclinations)
    [lat_s, si] = sort(results(k).lat_asc);
    dev_s = results(k).dev_asc(si);
    plot(lat_s, dev_s, '-', 'Color', colors{k}, 'LineWidth', 1.8, 'DisplayName', labels{k});
end
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('Geodetic Latitude (°)');
ylabel('Geodetic alt − spherical alt (km)');
title('WGS84 deviation – ascending pass');
legend('Location', 'best');
grid on;

%% ── Save ─────────────────────────────────────────────────────────────────
script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

exportgraphics(f1, fullfile(out_dir, 'altitude_vs_latitude_full_orbit.png'),    'Resolution', 300);
exportgraphics(f2, fullfile(out_dir, 'altitude_vs_latitude_ascending_pass.png'), 'Resolution', 300);
fprintf('\nSaved to %s\n', out_dir);
end
