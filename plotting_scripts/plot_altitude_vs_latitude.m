function plot_altitude_vs_latitude(HEIGHT_KM)
% PLOT_ALTITUDE_VS_LATITUDE
%   Propagates a single polar satellite with 'numerical' (J2+J3+J4+drag+SRP)
%   and 'sgp4', showing geodetic altitude (= nadir range) vs latitude and
%   vs time for one full orbit.
%
%   The sawtooth artefact of sorting by latitude is avoided by splitting
%   the trace into ascending (dlat>0) and descending (dlat<0) passes.
%
%   Usage:
%     plot_altitude_vs_latitude()       % default: 1000 km
%     plot_altitude_vs_latitude(600)

if nargin < 1 || isempty(HEIGHT_KM), HEIGHT_KM = 1000; end


%% ── Propagate ────────────────────────────────────────────────────────────
r_earth = 6378.14e3;              % matches coverage_simulator_function.m
a       = r_earth + HEIGHT_KM*1e3;
T_s     = 2*pi * sqrt(a^3 / 3.986004418e14);

sc            = satelliteScenario;
sc.StartTime  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
sc.StopTime   = sc.StartTime + seconds(T_s*4 + 5);
sc.SampleTime = 5;

sat_num = satellite(sc, a, 0, 90, 0, 0, 0, 'OrbitPropagator','numerical', 'Name','Numerical');
sat_sgp = satellite(sc, a, 0, 90, 0, 0, 0, 'OrbitPropagator','sgp4',     'Name','SGP4');

fprintf('Propagating  (H=%d km, T=%.1f min, dt=%d s)...', HEIGHT_KM, T_s/60, sc.SampleTime);
tic;
pos_num = squeeze(states(sat_num, 'CoordinateFrame','ECEF'))';   % [nT×3] m
pos_sgp = squeeze(states(sat_sgp, 'CoordinateFrame','ECEF'))';
fprintf(' %.1f s\n', toc);

%% ── Geodetic altitude via ecef2lla ───────────────────────────────────────
lla_num = ecef2lla(pos_num);  lat_num = lla_num(:,1);  alt_num = lla_num(:,3)/1e3;
lla_sgp = ecef2lla(pos_sgp);  lat_sgp = lla_sgp(:,1);  alt_sgp = lla_sgp(:,3)/1e3;

nT  = size(pos_num, 1);
t_m = (0:nT-1)' * sc.SampleTime / 60;   % elapsed time [minutes]

% Split ascending / descending to avoid sawtooth when plotting vs latitude
asc_n = [false; diff(lat_num) >= 0];
asc_s = [false; diff(lat_sgp) >= 0];

%% ── Summary ──────────────────────────────────────────────────────────────
n_qtr = round((T_s/4) / sc.SampleTime);
fprintf('\n%-12s  Equator  N-pole  S-pole  Range\n', 'Propagator');
fprintf('%-12s  %7.2f  %7.2f  %7.2f  %5.2f km\n', 'Numerical', ...
    alt_num(1), alt_num(n_qtr), alt_num(3*n_qtr), range(alt_num));
fprintf('%-12s  %7.2f  %7.2f  %7.2f  %5.2f km\n', 'SGP4', ...
    alt_sgp(1), alt_sgp(n_qtr), alt_sgp(3*n_qtr), range(alt_sgp));

%% ── Plot ─────────────────────────────────────────────────────────────────
set(0,'DefaultAxesFontSize',13,'DefaultTextFontSize',13);
c_num = '#0072BD';
c_sgp = '#D95319';

f = figure('Name', sprintf('Propagator Comparison  –  %d km', HEIGHT_KM), ...
    'Color','w', 'Position',[60 80 1200 460]);

%% Left: altitude vs latitude — ascending (solid) and descending (dashed)
subplot(1,2,1); hold on;
plot(lat_num(asc_n),  alt_num(asc_n),  '-',  'Color',c_num, 'LineWidth',2.0, 'DisplayName','Numerical asc');
plot(lat_num(~asc_n), alt_num(~asc_n), '--', 'Color',c_num, 'LineWidth',1.3, 'DisplayName','Numerical desc');
plot(lat_sgp(asc_s),  alt_sgp(asc_s),  '-',  'Color',c_sgp, 'LineWidth',2.0, 'DisplayName','SGP4 asc');
plot(lat_sgp(~asc_s), alt_sgp(~asc_s), '--', 'Color',c_sgp, 'LineWidth',1.3, 'DisplayName','SGP4 desc');
yline(HEIGHT_KM, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2, ...
    'DisplayName', sprintf('Nominal %d km', HEIGHT_KM));
xlabel('Geodetic latitude (°)');
ylabel('Geodetic altitude / nadir range (km)');
title('Altitude vs latitude', 'FontWeight','bold');
legend('Location','north', 'NumColumns',2);
grid on; xlim([-90 90]);

%% Right: altitude vs time
subplot(1,2,2); hold on;
plot(t_m, alt_num, '-',  'Color',c_num, 'LineWidth',2.0, 'DisplayName','Numerical');
plot(t_m, alt_sgp, '--', 'Color',c_sgp, 'LineWidth',1.6, 'DisplayName','SGP4');
yline(HEIGHT_KM, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2, ...
    'DisplayName', sprintf('Nominal %d km', HEIGHT_KM));
xlabel('Time (min)');
ylabel('Geodetic altitude / nadir range (km)');
title('Altitude vs time  (one orbit)', 'FontWeight','bold');
legend('Location','best');
grid on; xlim([0 t_m(end)]);

sgtitle(sprintf('Numerical vs SGP4  —  %d km orbit', HEIGHT_KM), ...
    'FontWeight','bold', 'FontSize',14);

%% Save
script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, 'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
fname = fullfile(out_dir, sprintf('propagator_comparison_%dkm.png', HEIGHT_KM));
exportgraphics(f, fname, 'Resolution',300);
fprintf('\nSaved → %s\n', fname);
end
