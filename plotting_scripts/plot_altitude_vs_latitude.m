function plot_altitude_vs_latitude(HEIGHT_KM_vec, el_min_vec, inc_deg)
% PLOT_ALTITUDE_VS_LATITUDE
%   Shows what coverage_simulator_function.m actually experiences when you
%   specify "H km circular orbit, e=0".
%
%   The precision simulator places satellites at:
%       r_sat = 6378.14 km  +  H km      (spherical orbit, constant radius)
%   and places UEs on the WGS84 ellipsoid via lla2ecef(alt=0).
%   Elevation angles are computed relative to the WGS84 surface normal.
%
%   Two effects shift the effective nadir slant range away from H:
%
%   (1) WGS84 GEOMETRY (~21 km, dominant, always present)
%       The WGS84 polar radius is 21.385 km shorter than the equatorial
%       radius.  The satellite orbit is spherical; the Earth surface is not.
%       Nadir range increases from 0 km extra at equator to +21.4 km at poles.
%
%   (2) J2 ORBITAL PERTURBATION (~8-9 km, correlated with latitude)
%       Earth's equatorial bulge deflects the circular orbit: the satellite
%       flies ~9 km LOWER at the equator and ~9 km HIGHER at the poles.
%       Analytical: dr = -(3/2) J2 (R_E/a)^2 a sin^2(i) cos(2*lat)
%       Compounds with WGS84 at poles; partially cancels at equator.
%
%   Combined:  equator  ~  H - 9 km   (J2 reduces nadir range)
%              poles    ~  H + 30 km  (WGS84 +21 km, J2 +9 km)
%
%   Usage:
%     plot_altitude_vs_latitude()
%     plot_altitude_vs_latitude([600 1000 1200], [10 25], 87)

if nargin < 1 || isempty(HEIGHT_KM_vec), HEIGHT_KM_vec = [600, 1000, 1200]; end
if nargin < 2 || isempty(el_min_vec),    el_min_vec    = [10, 25];          end
if nargin < 3 || isempty(inc_deg),       inc_deg       = 87;                end

%% Constants matching coverage_simulator_function.m exactly
r_earth = 6378.14e3;    % r_earth used in simulator  [m]
b_wgs   = 6356.752e3;   % WGS84 polar semi-minor axis [m]
J2      = 1.08263e-3;

fprintf('r_earth = %.3f km  (as in coverage_simulator_function.m)\n', r_earth/1e3);
fprintf('WGS84 equatorial - polar = %.3f km\n', (r_earth - b_wgs)/1e3);

%% Latitude grid
lat_deg = linspace(-90, 90, 3601);
phi     = deg2rad(lat_deg);

% Geocentric radius of WGS84 ellipsoid at geodetic latitude
r_wgs84 = sqrt( ((r_earth^2 .* cos(phi)).^2 + (b_wgs^2 .* sin(phi)).^2) ./ ...
                ((r_earth   .* cos(phi)).^2 + (b_wgs   .* sin(phi)).^2) );

eq_idx   = ceil(numel(lat_deg)/2);
pole_idx = numel(lat_deg);

%% Figure 1: Nadir slant range vs latitude
set(0,'DefaultAxesFontSize',13,'DefaultTextFontSize',13);
f1 = figure('Name','Precision Simulator: Effective Nadir Range', ...
    'Color','w','Position',[60 80 1100 480]);
hold on;
colors = lines(numel(HEIGHT_KM_vec));

fprintf('\n%-10s  %-14s  %-14s  %-14s  %-14s\n', ...
    'Height', 'Equator(WGS)', 'Pole(WGS)', 'Equator(+J2)', 'Pole(+J2)');

for ki = 1:numel(HEIGHT_KM_vec)
    H_m   = HEIGHT_KM_vec(ki) * 1e3;
    r_sat = r_earth + H_m;

    % WGS84 geometry: nadir range = r_sat - r_WGS84(lat)
    nadir_wgs = (r_sat - r_wgs84) / 1e3;

    % J2 analytical short-period radial deviation (approximate: lat ~ arg-of-lat)
    % dr = -(3/2) J2 (R_E/a)^2 a sin^2(i) cos(2*lat)
    dr_j2    = -(3/2) * J2 * (r_earth/r_sat)^2 * r_sat * sind(inc_deg)^2 .* cos(2*phi);
    nadir_j2 = nadir_wgs + dr_j2/1e3;

    fprintf('%-10s  %-14.1f  %-14.1f  %-14.1f  %-14.1f\n', ...
        sprintf('%d km', HEIGHT_KM_vec(ki)), ...
        nadir_wgs(eq_idx), nadir_wgs(pole_idx), ...
        nadir_j2(eq_idx),  nadir_j2(pole_idx));

    % Solid = WGS84 geometry only
    plot(lat_deg, nadir_wgs, '-', 'Color', colors(ki,:), 'LineWidth', 2.2, ...
        'DisplayName', sprintf('H=%d km  (WGS84 geom)', HEIGHT_KM_vec(ki)));

    % Dashed = WGS84 + J2
    plot(lat_deg, nadir_j2, '--', 'Color', colors(ki,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('H=%d km  (+J2, i=%d\xB0)', HEIGHT_KM_vec(ki), inc_deg));

    % Shaded band between the two
    fill([lat_deg, fliplr(lat_deg)], [nadir_wgs, fliplr(nadir_j2)], ...
        colors(ki,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Nominal horizontal reference
    yline(HEIGHT_KM_vec(ki), ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
end

text(-83, HEIGHT_KM_vec(1) - 12, 'Nominal (specified)', ...
    'Color', [0.55 0.55 0.55], 'FontSize', 10);

xlabel('Geodetic latitude (\circ)');
ylabel('Nadir slant range (km)');
title(sprintf(['Effective nadir range in coverage\\_simulator\\_function.m\n' ...
    'Solid = WGS84 geometry only   |   Dashed = +J2 (i = %d\xB0)'], inc_deg), ...
    'FontWeight','bold');
legend('Location','north','NumColumns',2);
grid on;

%% Figure 2: Coverage footprint half-angle vs latitude
%   rho = acos(r_WGS84(lat) * cos(el_min) / r_sat) - el_min
H_ref  = HEIGHT_KM_vec(ceil(numel(HEIGHT_KM_vec)/2));
r_ref  = r_earth + H_ref * 1e3;
colors2 = lines(numel(el_min_vec));

f2 = figure('Name', sprintf('Footprint Half-Angle  (H=%d km)', H_ref), ...
    'Color','w','Position',[60 580 900 420]);
hold on;

fprintf('\nCoverage half-angle  (H=%d km, WGS84 geometry):\n', H_ref);
fprintf('%-14s  %-12s  %-12s  %-12s\n', 'Min elev', 'Equator', 'Pole', 'Pole-Equator');

for ki = 1:numel(el_min_vec)
    eps = el_min_vec(ki);
    rho = acosd(min(1, r_wgs84 .* cosd(eps) / r_ref)) - eps;

    fprintf('%-14s  %-12.3f  %-12.3f  %+.3f\n', ...
        sprintf('%d\xB0', eps), rho(eq_idx), rho(pole_idx), rho(pole_idx)-rho(eq_idx));

    plot(lat_deg, rho, '-', 'Color', colors2(ki,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Min elev = %d\xB0', eps));
end

xlabel('Geodetic latitude (\circ)');
ylabel('Coverage half-angle \rho (\circ)');
title(sprintf('Coverage footprint half-angle vs latitude  (H = %d km)', H_ref), ...
    'FontWeight','bold');
legend('Location','north');
grid on;

%% Save
script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
exportgraphics(f1, fullfile(out_dir, 'nadir_range_vs_latitude.png'),       'Resolution',300);
exportgraphics(f2, fullfile(out_dir, 'coverage_halfangle_vs_latitude.png'), 'Resolution',300);
fprintf('\nSaved to %s\n', out_dir);
end
