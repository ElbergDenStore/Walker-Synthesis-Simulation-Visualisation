function plot_altitude_vs_latitude(HEIGHT_KM, el_min_vec)
% PLOT_ALTITUDE_VS_LATITUDE
%   Propagates a single polar satellite with SGP4 — the same pipeline as
%   coverage_simulator_function.m — and shows the resulting geodetic
%   altitude (= nadir slant range to WGS84 surface) vs latitude.
%
%   Uses exactly the same r_earth, propagator, and ecef2lla call as the
%   precision simulator.  A 90° polar orbit is used so the full latitude
%   range (-90° → +90°) is visible in one orbit.
%
%   Usage:
%     plot_altitude_vs_latitude()           % 1000 km, el=[10 25]
%     plot_altitude_vs_latitude(600)
%     plot_altitude_vs_latitude(1000, [5 10 20 35])

if nargin < 1 || isempty(HEIGHT_KM),  HEIGHT_KM  = 1000;     end
if nargin < 2 || isempty(el_min_vec), el_min_vec = [10, 25]; end

%% ── Propagate with SGP4 (identical to coverage_simulator_function.m) ────
r_earth = 6378.14e3;          % matches simulator
a       = r_earth + HEIGHT_KM*1e3;
T_s     = 2*pi * sqrt(a^3 / 3.986004418e14);   % approx period for stop time

sc            = satelliteScenario;
sc.StartTime  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
sc.StopTime   = sc.StartTime + seconds(T_s + 5);
sc.SampleTime = 5;

% Single satellite: 90° polar, e=0, nu=0 at equator northbound
sat = satellite(sc, a, 0, 90, 0, 0, 0, ...
    'OrbitPropagator','sgp4','Name','PolarSat');

fprintf('Propagating SGP4  (H=%d km, T=%.0f s, dt=%d s)...\n', HEIGHT_KM, T_s, 5);
[pos_raw, ~, ~] = states(sat, 'CoordinateFrame','ECEF');
pos = squeeze(pos_raw)';     % [nT x 3] metres

%% ── Geodetic altitude via ecef2lla (same call as simulator) ─────────────
lla      = ecef2lla(pos);
lat      = lla(:,1);          % geodetic latitude [deg]
alt_km   = lla(:,3) / 1e3;   % geodetic altitude = nadir slant range [km]

% Ascending pass only (equator northbound → north pole): first T/4
n_asc = round((T_s/4) / sc.SampleTime);
asc   = 1:n_asc;

fprintf('\nSGP4 nadir slant range (geodetic altitude from ecef2lla):\n');
fprintf('  Equator (ascending) : %.2f km\n', alt_km(1));
fprintf('  North pole          : %.2f km\n', alt_km(n_asc));
fprintf('  South pole          : %.2f km\n', alt_km(3*n_asc));
fprintf('  Pole - equator      : %+.2f km\n', alt_km(n_asc) - alt_km(1));
fprintf('  Peak-to-peak        : %.2f km\n', range(alt_km));

%% ── Coverage half-angle from actual SGP4 radius ─────────────────────────
%   rho = acos(r_WGS84(lat) * cos(el) / r_sat)  -- but we use the real r_sat
%   from SGP4 positions rather than nominal a.
r_sat_vec = sqrt(sum(pos.^2, 2));   % actual SGP4 orbital radius [m]

b_wgs   = 6356.752e3;  % WGS84 polar radius
phi     = deg2rad(lat);
r_wgs84 = sqrt( ((r_earth^2.*cos(phi)).^2 + (b_wgs^2.*sin(phi)).^2) ./ ...
                ((r_earth  .*cos(phi)).^2 + (b_wgs  .*sin(phi)).^2) );

rho = zeros(numel(el_min_vec), numel(lat));
for ki = 1:numel(el_min_vec)
    arg = r_wgs84 .* cosd(el_min_vec(ki)) ./ r_sat_vec;
    rho(ki,:) = acosd(min(1, arg)) - el_min_vec(ki);
end

%% ── Plot ─────────────────────────────────────────────────────────────────
set(0,'DefaultAxesFontSize',13,'DefaultTextFontSize',13);
colors = lines(numel(el_min_vec));

f = figure('Name', sprintf('SGP4 Simulator Geometry  –  %d km orbit', HEIGHT_KM), ...
    'Color','w','Position',[60 80 1200 460]);

%% Left: geodetic altitude (nadir slant range) vs latitude — full orbit
subplot(1,2,1);
[lat_s, si] = sort(lat);
plot(lat_s, alt_km(si), '-', 'Color','#0072BD', 'LineWidth',2.2);
hold on;
yline(HEIGHT_KM, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1.4, ...
    'DisplayName', sprintf('Nominal %d km', HEIGHT_KM));

% Annotate pole deviation
pole_dev = alt_km(n_asc) - alt_km(1);
text(60, alt_km(n_asc) - pole_dev*0.35, ...
    sprintf('+%.1f km\nat poles', pole_dev), ...
    'FontSize',11,'Color','#0072BD','HorizontalAlignment','center');

xlabel('Geodetic latitude (°)');
ylabel('Geodetic altitude / nadir range (km)');
title(sprintf('SGP4 nadir slant range  (H = %d km nominal)', HEIGHT_KM), ...
    'FontWeight','bold');
legend({'SGP4 + ecef2lla', sprintf('Nominal %d km', HEIGHT_KM)}, 'Location','north');
grid on; xlim([-90 90]);

%% Right: coverage half-angle vs latitude — ascending pass
subplot(1,2,2);
hold on;
for ki = 1:numel(el_min_vec)
    plot(lat(asc), rho(ki,asc), '-', 'Color',colors(ki,:), 'LineWidth',2, ...
        'DisplayName', sprintf('Min el = %d°', el_min_vec(ki)));
end
xlabel('Geodetic latitude (°)');
ylabel('Coverage half-angle \rho (°)');
title('Footprint half-angle  (ascending pass)', 'FontWeight','bold');
legend('Location','southeast');
grid on; xlim([0 90]);

sgtitle(sprintf('coverage\\_simulator\\_function.m  —  %d km orbit, SGP4 propagation', HEIGHT_KM), ...
    'FontWeight','bold','FontSize',14);

%% Save
script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, 'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
fname = fullfile(out_dir, sprintf('simulator_geometry_%dkm.png', HEIGHT_KM));
exportgraphics(f, fname, 'Resolution',300);
fprintf('\nSaved to %s\n', fname);
end
