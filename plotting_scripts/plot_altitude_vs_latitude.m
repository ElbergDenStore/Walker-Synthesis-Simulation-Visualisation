function plot_altitude_vs_latitude()
% PLOT_ALTITUDE_VS_LATITUDE
%   Illustrates the WGS84 geometric altitude effect for a 90° polar orbit
%   at 1000 km nominal height.
%
%   WHY IS ALTITUDE HIGHEST AT THE POLES?
%   The WGS84 ellipsoid has an equatorial radius of 6378.137 km but a polar
%   radius of only 6356.752 km — a difference of ~21.4 km.  A satellite
%   flying at constant orbital radius therefore sits ~21.4 km farther above
%   the WGS84 surface at the poles than at the equator.  This effect is
%   perfectly symmetric about the equator (same deviation at ±90°).
%
%   The WGS84 deviation is isolated by locking |r| = a_design (removing
%   any real orbit perturbation) and converting with ecef2lla.
%
%   The "combined" curve uses real SGP4 positions with ecef2lla.  It may
%   show a small north-south asymmetry (~1-3 km) from J2 perturbation, but
%   the dominant ~21 km effect is purely geometric.
%
%   Figure 1: full orbital period    (all latitudes)
%   Figure 2: ascending pass only    (equator → north pole)

HEIGHT_KM = 1000;
r_sphere  = 6378.14e3;          % spherical Earth used for orbit design [m]
a_design  = r_sphere + HEIGHT_KM*1e3;

%% ── 1. Propagate one full orbit with SGP4 ───────────────────────────────
T_orbit_s  = 2*pi * sqrt(a_design^3 / 3.986004418e14);
sample_s   = 5;

sc            = satelliteScenario;
sc.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
sc.StopTime   = sc.StartTime + seconds(T_orbit_s + sample_s);
sc.SampleTime = sample_s;

sat = satellite(sc, a_design, 0, 90, 0, 0, 0, ...
    'OrbitPropagator', 'sgp4', 'Name', 'PolarSat');

fprintf('Propagating (%.0f s period, %d s sample)...\n', T_orbit_s, sample_s);
[pos_raw, ~, ~] = states(sat, "CoordinateFrame", "ECEF");
pos = squeeze(pos_raw)';   % [nT x 3] metres

%% ── 2. Compute altitude effects ─────────────────────────────────────────
r_actual = sqrt(sum(pos.^2, 2));   % real SGP4 orbital radius [m]

% COMBINED: real SGP4 orbit + WGS84 geodetic conversion
lla_real     = ecef2lla(pos);
lat          = lla_real(:,1);
alt_combined = lla_real(:,3) / 1e3;

% WGS84 ONLY: lock |r| = a_design (no J2 perturbation), then apply ecef2lla.
% This isolates the purely geometric effect of the oblate reference ellipsoid:
% the WGS84 polar radius (6356.75 km) is ~21.4 km shorter than the equatorial
% radius, so a satellite at constant orbital radius reads higher altitude over
% the poles.  Effect is perfectly symmetric about the equator.
pos_const_r    = pos ./ r_actual .* a_design;
lla_const      = ecef2lla(pos_const_r);
alt_wgs84_only = lla_const(:,3) / 1e3;

% Deviations from nominal HEIGHT_KM (always >= 0: WGS84 surface dips below
% the sphere at all latitudes except the equator)
dev_wgs84    = alt_wgs84_only - HEIGHT_KM;
dev_combined = alt_combined   - HEIGHT_KM;

% Ascending pass: nu=0 at epoch → equator northbound.  First T/4 = pole.
n_asc = round((T_orbit_s/4) / sample_s);
asc   = 1:n_asc;

fprintf('WGS84    p-p: %.4f km  (symmetric; WGS84 polar radius ~21.4 km below equatorial)\n', range(dev_wgs84));
fprintf('Combined p-p: %.4f km  (small asymmetry from real J2 orbital perturbation)\n', range(dev_combined));

%% ── 3. Plotting ──────────────────────────────────────────────────────────
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontSize', 13);
c_wgs  = '#0072BD';
c_comb = '#000000';

    function draw_panel(lat_v, devs, names, clrs, ttl)
        hold on;
        for ki = 1:numel(devs)
            [ls, si] = sort(lat_v);
            plot(ls, devs{ki}(si), '-', 'Color', clrs{ki}, ...
                'LineWidth', 1.8, 'DisplayName', names{ki});
        end
        yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
        xlabel('Geodetic latitude (°)');
        ylabel('Deviation from nominal (km)');
        title(ttl);
        legend('Location', 'best');
        grid on;
    end

%% Figure 1 – full orbit
f1 = figure('Name', 'Altitude Effects – Full Orbit', ...
    'Color', 'w', 'Position', [60 80 1100 430]);

subplot(1,2,1);
draw_panel(lat, {dev_wgs84}, {'WGS84 geometric'}, {c_wgs}, ...
    sprintf('WGS84 effect  (p-p: %.3f km)', range(dev_wgs84)));

subplot(1,2,2);
draw_panel(lat, {dev_wgs84, dev_combined}, ...
    {'WGS84 geometric', 'Combined (SGP4 + WGS84)'}, ...
    {c_wgs, c_comb}, 'Combined effect');

sgtitle(sprintf('Altitude variation – 90° polar orbit, %d km nominal (full orbit)', HEIGHT_KM), ...
    'FontWeight', 'bold', 'FontSize', 14);

%% Figure 2 – ascending pass
f2 = figure('Name', 'Altitude Effects – Ascending Pass', ...
    'Color', 'w', 'Position', [60 560 1100 430]);

subplot(1,2,1);
draw_panel(lat(asc), {dev_wgs84(asc)}, {'WGS84 geometric'}, {c_wgs}, ...
    sprintf('WGS84 effect  (p-p: %.3f km)', range(dev_wgs84(asc))));

subplot(1,2,2);
draw_panel(lat(asc), {dev_wgs84(asc), dev_combined(asc)}, ...
    {'WGS84 geometric', 'Combined (SGP4 + WGS84)'}, ...
    {c_wgs, c_comb}, 'Combined effect');

sgtitle(sprintf('Altitude variation – 90° polar orbit, %d km nominal (ascending pass)', HEIGHT_KM), ...
    'FontWeight', 'bold', 'FontSize', 14);

%% ── Save ─────────────────────────────────────────────────────────────────
script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

exportgraphics(f1, fullfile(out_dir, 'altitude_vs_latitude_full_orbit.png'),     'Resolution', 300);
exportgraphics(f2, fullfile(out_dir, 'altitude_vs_latitude_ascending_pass.png'),  'Resolution', 300);
fprintf('Saved to %s\n', out_dir);
end
