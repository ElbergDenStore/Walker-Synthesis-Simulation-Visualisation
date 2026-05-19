function plot_altitude_vs_latitude()
% PLOT_ALTITUDE_VS_LATITUDE
%   Uses the MATLAB Satellite Toolbox (SGP4) to separate and illustrate the
%   two sources of altitude variation for a circular 90° polar orbit at
%   1000 km nominal height:
%
%     (1) WGS84 geometric effect  – oblate reference ellipsoid makes a
%         satellite at constant orbital radius read ~21 km higher at poles.
%         Isolated by locking |r| = a_design and converting with ecef2lla.
%
%     (2) J2 orbital perturbation – Earth's equatorial bulge perturbs the
%         circular orbit, causing real ~1-3 km variation in orbital radius.
%         Isolated by using the real SGP4 |r| against a perfect sphere.
%
%     (3) Combined – real SGP4 positions converted with ecef2lla.
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

%% ── 2. Separate the two effects ─────────────────────────────────────────
r_actual = sqrt(sum(pos.^2, 2));   % real SGP4 orbital radius [m]

% COMBINED: real SGP4 orbit + WGS84 geodetic conversion
lla_real     = ecef2lla(pos);
lat          = lla_real(:,1);
alt_combined = lla_real(:,3) / 1e3;

% WGS84 ONLY: lock |r| = a_design (no J2), then apply ecef2lla.
% Shows purely how the oblate reference ellipsoid changes the altitude reading
% for a satellite at constant orbital radius.
pos_const_r    = pos ./ r_actual .* a_design;
lla_const      = ecef2lla(pos_const_r);
alt_wgs84_only = lla_const(:,3) / 1e3;

% J2 ONLY: deviation of the actual orbital radius from its orbit-mean.
% Subtracting mean(r_actual) instead of a_design removes the DC offset
% introduced by SGP4's osculating-to-mean element conversion at epoch
% (a_mean ≈ a_design + ~4 km), leaving only the real short-period J2
% oscillation centred on zero.
r_mean      = mean(r_actual);
alt_j2_only = (r_actual - r_mean) / 1e3;   % km, zero-mean oscillation

% Reference deviations
% Note: dev_wgs84 and dev_combined are measured from the spherical design
% altitude. Both are always >= 0 because the WGS84 surface dips below the
% sphere (sphere radius = WGS84 equatorial radius) at all non-zero latitudes.
dev_wgs84    = alt_wgs84_only - HEIGHT_KM;
dev_j2       = alt_j2_only;                 % already zero-mean
dev_combined = alt_combined  - HEIGHT_KM;

% Ascending pass: nu=0 at epoch → equator northbound.  First T/4 = pole.
n_asc = round((T_orbit_s/4) / sample_s);
asc   = 1:n_asc;

fprintf('WGS84  p-p: %.4f km  (always >= 0: WGS84 surface below sphere)\n', range(dev_wgs84));
fprintf('J2     p-p: %.4f km  (zero-mean oscillation)\n', range(dev_j2));
fprintf('J2 mean offset from a_design: %.3f km  (SGP4 osc->mean conversion)\n', (r_mean - a_design)/1e3);
fprintf('Combined p-p: %.4f km\n', range(dev_combined));

%% ── 3. Plotting ──────────────────────────────────────────────────────────
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontSize', 13);
c_wgs  = '#0072BD';
c_j2   = '#D95319';
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
    'Color', 'w', 'Position', [60 80 1500 430]);

subplot(1,3,1);
draw_panel(lat, {dev_wgs84}, {'WGS84 geometric'}, {c_wgs}, ...
    sprintf('WGS84 effect  (p-p: %.3f km)', range(dev_wgs84)));

subplot(1,3,2);
draw_panel(lat, {dev_j2}, {'J2 perturbation'}, {c_j2}, ...
    sprintf('J2 perturbation  (p-p: %.4f km)', range(dev_j2)));

subplot(1,3,3);
draw_panel(lat, {dev_wgs84, dev_j2, dev_combined}, ...
    {'WGS84 geometric', 'J2 perturbation', 'Combined'}, ...
    {c_wgs, c_j2, c_comb}, 'All effects combined');

sgtitle(sprintf('Altitude variation – 90° polar orbit, %d km nominal (full orbit)', HEIGHT_KM), ...
    'FontWeight', 'bold', 'FontSize', 14);

%% Figure 2 – ascending pass
f2 = figure('Name', 'Altitude Effects – Ascending Pass', ...
    'Color', 'w', 'Position', [60 560 1500 430]);

subplot(1,3,1);
draw_panel(lat(asc), {dev_wgs84(asc)}, {'WGS84 geometric'}, {c_wgs}, ...
    sprintf('WGS84 effect  (p-p: %.3f km)', range(dev_wgs84(asc))));

subplot(1,3,2);
draw_panel(lat(asc), {dev_j2(asc)}, {'J2 perturbation'}, {c_j2}, ...
    sprintf('J2 perturbation  (p-p: %.4f km)', range(dev_j2(asc))));

subplot(1,3,3);
draw_panel(lat(asc), {dev_wgs84(asc), dev_j2(asc), dev_combined(asc)}, ...
    {'WGS84 geometric', 'J2 perturbation', 'Combined'}, ...
    {c_wgs, c_j2, c_comb}, 'All effects combined');

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
