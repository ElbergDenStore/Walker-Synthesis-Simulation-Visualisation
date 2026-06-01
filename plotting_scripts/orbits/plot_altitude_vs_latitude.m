function plot_altitude_vs_latitude(HEIGHT_KM)
% PLOT_ALTITUDE_VS_LATITUDE
%   Propagates a single polar satellite with 'two-body-keplerian', 'numerical'
%   (J2+J3+J4+drag+SRP), and 'sgp4', showing geodetic altitude (= nadir range)
%   vs latitude and vs time for one full orbit.
%
%   The sawtooth artefact of sorting by latitude is avoided by splitting
%   the trace into ascending (dlat>0) and descending (dlat<0) passes.
%
%   Usage:
%     plot_altitude_vs_latitude()       % default: 1000 km
%     plot_altitude_vs_latitude(600)

if nargin < 1 || isempty(HEIGHT_KM), HEIGHT_KM = 1000; end


%% ── Propagate ────────────────────────────────────────────────────────────
r_earth = 6378.14e3;              % matches constellation_simulator.m
a       = r_earth + HEIGHT_KM*1e3;
T_s     = 2*pi * sqrt(a^3 / 3.986004418e14);

sc            = satelliteScenario;
sc.StartTime  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
sc.StopTime   = sc.StartTime + seconds(T_s*2 + 5);
sc.SampleTime = 5;

sat_kepler = satellite(sc, a, 0, 90, 0, 0, 0, 'OrbitPropagator','two-body-keplerian', 'Name','Kepler');
sat_num    = satellite(sc, a, 0, 90, 0, 0, 0, 'OrbitPropagator','numerical',          'Name','Numerical');
sat_sgp    = satellite(sc, a, 0, 90, 0, 0, 0, 'OrbitPropagator','sgp4',               'Name','SGP4');

fprintf('Propagating  (H=%d km, T=%.1f min, dt=%d s)...', HEIGHT_KM, T_s/60, sc.SampleTime);
tic;
pos_kep = squeeze(states(sat_kepler, 'CoordinateFrame','ECEF'))';  % [nT×3] m
pos_num = squeeze(states(sat_num,    'CoordinateFrame','ECEF'))';
pos_sgp = squeeze(states(sat_sgp,    'CoordinateFrame','ECEF'))';
fprintf(' %.1f s\n', toc);

%% ── Geodetic altitude via ecef2lla ───────────────────────────────────────
lla_kep = ecef2lla(pos_kep);  lat_kep = lla_kep(:,1);  alt_kep = lla_kep(:,3)/1e3;
lla_num = ecef2lla(pos_num);  lat_num = lla_num(:,1);  alt_num = lla_num(:,3)/1e3;
lla_sgp = ecef2lla(pos_sgp);  lat_sgp = lla_sgp(:,1);  alt_sgp = lla_sgp(:,3)/1e3;

nT  = size(pos_num, 1);
t_m = (0:nT-1)' * sc.SampleTime / 60;   % elapsed time [minutes]

%% ── Summary ──────────────────────────────────────────────────────────────
n_qtr = round((T_s/4) / sc.SampleTime);
fprintf('\n%-12s  Equator  N-pole  S-pole  Range\n', 'Propagator');
fprintf('%-12s  %7.2f  %7.2f  %7.2f  %5.2f km\n', 'Kepler', ...
    alt_kep(1), alt_kep(n_qtr), alt_kep(3*n_qtr), range(alt_kep));
fprintf('%-12s  %7.2f  %7.2f  %7.2f  %5.2f km\n', 'Numerical', ...
    alt_num(1), alt_num(n_qtr), alt_num(3*n_qtr), range(alt_num));
fprintf('%-12s  %7.2f  %7.2f  %7.2f  %5.2f km\n', 'SGP4', ...
    alt_sgp(1), alt_sgp(n_qtr), alt_sgp(3*n_qtr), range(alt_sgp));

%% ── Sort by latitude (removes sawtooth on lat axis) ─────────────────────
[lat_kep_s, si] = sort(lat_kep);  alt_kep_s = alt_kep(si);
[lat_num_s, si] = sort(lat_num);  alt_num_s = alt_num(si);
[lat_sgp_s, si] = sort(lat_sgp);  alt_sgp_s = alt_sgp(si);

%% ── Plot ─────────────────────────────────────────────────────────────────
set(0,'DefaultAxesFontSize',13,'DefaultTextFontSize',13);
c_kep = '#77AC30';
c_num = '#0072BD';
c_sgp = '#D95319';

script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, 'figures/propagators');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% Figure 1: altitude vs latitude (sorted) --------------------------------
f1 = figure('Name', sprintf('Altitude vs Latitude  –  %d km', HEIGHT_KM), ...
    'Color','w', 'Position',[60 80 500 400]);
hold on;
plot(lat_kep_s, alt_kep_s, '-',  'Color',c_kep, 'LineWidth',2.0, 'DisplayName','Kepler');
plot(lat_num_s, alt_num_s, '-',  'Color',c_num, 'LineWidth',2.0, 'DisplayName','Numerical');
plot(lat_sgp_s, alt_sgp_s, '-', 'Color',c_sgp, 'LineWidth',1.6, 'DisplayName','SGP4');
% yline(HEIGHT_KM, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2, ...
    % 'DisplayName', sprintf('Nominal %d km', HEIGHT_KM));
xlabel('Geodetic latitude (°)');
ylabel('Geodetic altitude / nadir range (km)');
title(sprintf('Altitude vs latitude  —  %d km orbit', HEIGHT_KM), 'FontWeight','bold');
legend('Location','north', 'NumColumns',2);
grid on; xlim([-90 90]);

fname1 = fullfile(out_dir, sprintf('propagator_alt_vs_lat_%dkm.png', HEIGHT_KM));
exportgraphics(f1, fname1, 'Resolution',300);
fprintf('\nSaved → %s\n', fname1);

%% Figure 2: altitude vs time ─────────────────────────────────────────────
f2 = figure('Name', sprintf('Altitude vs Time  –  %d km', HEIGHT_KM), ...
    'Color','w', 'Position',[100 80 500 400]);
hold on;
plot(t_m, alt_kep, '-',  'Color',c_kep, 'LineWidth',2.0, 'DisplayName','Kepler');
plot(t_m, alt_num, '-',  'Color',c_num, 'LineWidth',2.0, 'DisplayName','Numerical');
plot(t_m, alt_sgp, '-', 'Color',c_sgp, 'LineWidth',1.6, 'DisplayName','SGP4');
% yline(HEIGHT_KM, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2, ...
%     'DisplayName', sprintf('Nominal %d km', HEIGHT_KM));
xlabel('Time (min)');
ylabel('Geodetic altitude / nadir range (km)');
title(sprintf('Altitude vs time  —  %d km orbit', HEIGHT_KM), 'FontWeight','bold');
legend('Location','best');
grid on; xlim([0 t_m(end)]);

fname2 = fullfile(out_dir, sprintf('propagator_alt_vs_time_%dkm.png', HEIGHT_KM));
exportgraphics(f2, fname2, 'Resolution',300);
fprintf('Saved → %s\n', fname2);

%% Figure 3: deviation from nominal altitude vs time ──────────────────────
f3 = figure('Name', sprintf('Altitude Deviation  –  %d km', HEIGHT_KM), ...
    'Color','w', 'Position',[140 80 500 400]);
hold on;
plot(t_m, alt_kep - HEIGHT_KM, '-',  'Color',c_kep, 'LineWidth',2.0, 'DisplayName','Kepler');
plot(t_m, alt_num - HEIGHT_KM, '-',  'Color',c_num, 'LineWidth',2.0, 'DisplayName','Numerical');
plot(t_m, alt_sgp - HEIGHT_KM, '-', 'Color',c_sgp, 'LineWidth',1.6, 'DisplayName','SGP4');
% yline(0, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2, 'DisplayName','Nominal');
xlabel('Time (min)');
ylabel('Altitude deviation from nominal (km)');
title(sprintf('Deviation from nominal altitude  —  %d km orbit', HEIGHT_KM), 'FontWeight','bold');
legend('Location','best');
grid on; xlim([0 t_m(end)]);

fname3 = fullfile(out_dir, sprintf('propagator_alt_deviation_%dkm.png', HEIGHT_KM));
exportgraphics(f3, fname3, 'Resolution',300);
fprintf('Saved → %s\n', fname3);
end
