function plot_altitude_vs_latitude(HEIGHT_KM, el_min_vec)
% PLOT_ALTITUDE_VS_LATITUDE
%   Propagates a single polar satellite with both the 'numerical' propagator
%   (RK89, J2+J3+J4+drag+SRP) and 'sgp4', and overlays the results so the
%   SGP4 TLE-conversion artefact is visible against the physically rigorous
%   numerical solution.
%
%   r_earth and ecef2lla match coverage_simulator_function.m exactly.
%   A 90° polar orbit is used so the full latitude range is visible.
%
%   Usage:
%     plot_altitude_vs_latitude()           % 1000 km, el=[10 25]
%     plot_altitude_vs_latitude(600)
%     plot_altitude_vs_latitude(1000, [5 10 20 35])

if nargin < 1 || isempty(HEIGHT_KM),  HEIGHT_KM  = 1000;     end
if nargin < 2 || isempty(el_min_vec), el_min_vec = [10, 25]; end

%% ── Shared scenario parameters ───────────────────────────────────────────
r_earth = 6378.14e3;          % matches coverage_simulator_function.m
a       = r_earth + HEIGHT_KM*1e3;
T_s     = 2*pi * sqrt(a^3 / 3.986004418e14);   % Keplerian period [s]

sc            = satelliteScenario;
sc.StartTime  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
sc.StopTime   = sc.StartTime + seconds(T_s + 5);
sc.SampleTime = 5;

%% ── Propagate: numerical (J2+J3+J4+drag+SRP) ────────────────────────────
sat_num = satellite(sc, a, 0, 90, 0, 0, 0, ...
    'OrbitPropagator','numerical','Name','Numerical');

%% ── Propagate: SGP4 (same as coverage_simulator_function.m) ─────────────
sat_sgp = satellite(sc, a, 0, 90, 0, 0, 0, ...
    'OrbitPropagator','sgp4','Name','SGP4');

fprintf('Propagating both propagators  (H=%d km, T=%.0f s, dt=%d s)...\n', HEIGHT_KM, T_s, 5);
pos_num = squeeze(states(sat_num, 'CoordinateFrame','ECEF'))';   % [nT x 3] m
pos_sgp = squeeze(states(sat_sgp, 'CoordinateFrame','ECEF'))';

%% ── Geodetic altitude via ecef2lla ───────────────────────────────────────
lla_num  = ecef2lla(pos_num);
lat_num  = lla_num(:,1);
alt_num  = lla_num(:,3) / 1e3;

lla_sgp  = ecef2lla(pos_sgp);
lat_sgp  = lla_sgp(:,1);
alt_sgp  = lla_sgp(:,3) / 1e3;

n_asc = round((T_s/4) / sc.SampleTime);
asc   = 1:n_asc;

fprintf('\n%-30s  Equator    Pole    Pole-Eq\n', 'Propagator');
fprintf('%-30s  %7.2f  %7.2f  %+7.2f km\n', 'Numerical', ...
    alt_num(1), alt_num(n_asc), alt_num(n_asc)-alt_num(1));
fprintf('%-30s  %7.2f  %7.2f  %+7.2f km\n', 'SGP4 (simulator)', ...
    alt_sgp(1), alt_sgp(n_asc), alt_sgp(n_asc)-alt_sgp(1));

%% ── Coverage half-angle (numerical r_sat) ────────────────────────────────
r_sat_num = sqrt(sum(pos_num.^2, 2));
r_sat_sgp = sqrt(sum(pos_sgp.^2, 2));

b_wgs = 6356.752e3;
phi_n = deg2rad(lat_num);
r_wgs84_num = sqrt( ((r_earth^2.*cos(phi_n)).^2 + (b_wgs^2.*sin(phi_n)).^2) ./ ...
                    ((r_earth  .*cos(phi_n)).^2 + (b_wgs  .*sin(phi_n)).^2) );

phi_s = deg2rad(lat_sgp);
r_wgs84_sgp = sqrt( ((r_earth^2.*cos(phi_s)).^2 + (b_wgs^2.*sin(phi_s)).^2) ./ ...
                    ((r_earth  .*cos(phi_s)).^2 + (b_wgs  .*sin(phi_s)).^2) );

rho_num = zeros(numel(el_min_vec), numel(lat_num));
rho_sgp = zeros(numel(el_min_vec), numel(lat_sgp));
for ki = 1:numel(el_min_vec)
    rho_num(ki,:) = acosd(min(1, r_wgs84_num .* cosd(el_min_vec(ki)) ./ r_sat_num)) - el_min_vec(ki);
    rho_sgp(ki,:) = acosd(min(1, r_wgs84_sgp .* cosd(el_min_vec(ki)) ./ r_sat_sgp)) - el_min_vec(ki);
end

%% ── Plot ─────────────────────────────────────────────────────────────────
set(0,'DefaultAxesFontSize',13,'DefaultTextFontSize',13);
colors = lines(numel(el_min_vec));

f = figure('Name', sprintf('Propagator Comparison  –  %d km orbit', HEIGHT_KM), ...
    'Color','w','Position',[60 80 1200 460]);

%% Left: geodetic altitude vs latitude — full orbit, both propagators
subplot(1,2,1);
[lat_ns, si_n] = sort(lat_num);
[lat_ss, si_s] = sort(lat_sgp);
plot(lat_ns, alt_num(si_n), '-',  'Color','#0072BD', 'LineWidth',2.2, 'DisplayName','Numerical');
hold on;
plot(lat_ss, alt_sgp(si_s), '--', 'Color','#D95319', 'LineWidth',1.8, 'DisplayName','SGP4 (simulator)');
yline(HEIGHT_KM, ':', 'Color',[0.55 0.55 0.55], 'LineWidth',1.4, ...
    'DisplayName', sprintf('Nominal %d km', HEIGHT_KM));
xlabel('Geodetic latitude (°)');
ylabel('Geodetic altitude / nadir range (km)');
title(sprintf('Nadir slant range  (H = %d km nominal)', HEIGHT_KM), 'FontWeight','bold');
legend('Location','north');
grid on; xlim([-90 90]);

%% Right: coverage half-angle, ascending pass, numerical only
subplot(1,2,2);
hold on;
for ki = 1:numel(el_min_vec)
    plot(lat_num(asc), rho_num(ki,asc), '-',  'Color',colors(ki,:), 'LineWidth',2, ...
        'DisplayName', sprintf('Numerical  el=%d°', el_min_vec(ki)));
    plot(lat_sgp(asc), rho_sgp(ki,asc), '--', 'Color',colors(ki,:), 'LineWidth',1.4, ...
        'HandleVisibility','off');
end
% dummy entry for SGP4 line style in legend
plot(NaN, NaN, 'k--', 'LineWidth',1.4, 'DisplayName','SGP4 (dashed)');
xlabel('Geodetic latitude (°)');
ylabel('Coverage half-angle \rho (°)');
title('Footprint half-angle  (ascending pass)', 'FontWeight','bold');
legend('Location','southeast');
grid on; xlim([0 90]);

sgtitle(sprintf('Numerical vs SGP4  —  %d km orbit', HEIGHT_KM), ...
    'FontWeight','bold','FontSize',14);

%% Save
script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, 'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
fname = fullfile(out_dir, sprintf('simulator_geometry_%dkm_num_vs_sgp4.png', HEIGHT_KM));
exportgraphics(f, fname, 'Resolution',300);
fprintf('\nSaved to %s\n', fname);
end
