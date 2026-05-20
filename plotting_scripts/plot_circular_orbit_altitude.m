function plot_circular_orbit_altitude(HEIGHT_KM, INC_DEG)
% PLOT_CIRCULAR_ORBIT_ALTITUDE
%   Produces TWO publication-quality figures:
%
%   Figure 1 – Orbital altitude vs geodetic latitude (single pass).
%
%   Figure 2 – Same altitude trace on the primary y-axis, plus a secondary
%   y-axis with two contributing effects:
%     • WGS84 surface variation: how much the reference ellipsoid radius
%       varies with latitude (Earth's equatorial bulge).
%     • J2 orbit perturbation: how much the satellite's geocentric radius
%       deviates from the nominal semi-major axis due to J2.
%
%   Propagated with MATLAB's numerical propagator.
%
%   Usage:
%       plot_circular_orbit_altitude()              % defaults
%       plot_circular_orbit_altitude(1000, 90)

    if nargin < 1 || isempty(HEIGHT_KM), HEIGHT_KM = 1000; end
    if nargin < 2 || isempty(INC_DEG),   INC_DEG   = 90;   end

    %% Constants ------------------------------------------------------------
    Re_m   = 6378.137e3;   % WGS-84 equatorial (used for ECEF + WGS84 calc)
    b_WGS  = 6356752.3142; % WGS-84 polar radius (m)
    mu     = 3.986004418e14;
    J2     = 1.0826e-3;

    a_m  = Re_m + HEIGHT_KM*1e3;
    T_s  = 2*pi * sqrt(a_m^3 / mu);

    e_circ = (J2 / 2) * (Re_m / a_m)^2;
    w_deg  = 0;
    nu_deg = 0;

    fprintf('\nCircular-orbit altitude vs latitude\n');
    fprintf('  h = %g km   i = %g deg\n', HEIGHT_KM, INC_DEG);
    fprintf('  e_init = %.4e   Period: %.2f min\n\n', e_circ, T_s/60);

    %% Propagate ------------------------------------------------------------
    sc = satelliteScenario;
    sc.StartTime  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
    sc.StopTime   = sc.StartTime + seconds(T_s + 5);   % one orbit
    sc.SampleTime = 5;

    sat = satellite(sc, a_m, e_circ, INC_DEG, 0, w_deg, nu_deg, ...
        'Name', 'circular', 'OrbitPropagator', 'numerical');

    P   = squeeze(states(sat, 'CoordinateFrame','ECEF'))';  % Nx3, metres
    lla = ecef2lla(P);
    lat = lla(:,1);          % geodetic latitude (deg)
    alt = lla(:,3) / 1e3;   % geodetic altitude (km)

    % One clean pass: use ascending half, sorted by latitude
    asc_mask    = [false; diff(lat) >= 0];
    [lat_s, ix] = sort(lat(asc_mask));
    tmp         = alt(asc_mask);  alt_s = tmp(ix);

    %% Derived quantities for Figure 2 -------------------------------------
    % WGS84 geocentric surface radius at each latitude (analytical)
    lat_rad    = deg2rad(lat_s);
    R_wgs84_km = sqrt( ((Re_m^2 * cos(lat_rad)).^2 + (b_WGS^2 * sin(lat_rad)).^2) ./ ...
                       ((Re_m   * cos(lat_rad)).^2 + (b_WGS   * sin(lat_rad)).^2) ) / 1e3;
    wgs84_dev  = R_wgs84_km - mean(R_wgs84_km);   % deviation from its own mean (km)

    % Analytical first-order J2 short-period radial perturbation (zero-mean):
    %   δr = -(3/2) * a*(J2/2)*(Re/a)^2 * cos(2φ)
    % The raw formula (1 - 3sin²φ) has a DC offset of -1/2 that belongs to
    % the secular mean-motion correction, not the oscillatory perturbation.
    % Subtracting it (using cos(2φ) = 1 - 2sin²φ) gives the symmetric result:
    %   negative near the equator (J2 pulls orbit inward),
    %   positive near the poles  (J2 pushes orbit outward).
    j2_dev = -(3/2) * (a_m/1e3) * (J2/2) * (Re_m/a_m)^2 * cos(2*lat_rad);

    %% ---- Figure 1: simple altitude trace --------------------------------
    script_dir = fileparts(mfilename('fullpath'));
    out_dir    = fullfile(script_dir, 'figures');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    f1 = figure('Color','w', 'Position',[80 120 760 420]);
    hold on;
    plot(lat_s, alt_s, '-', 'Color',[0.00 0.45 0.74], 'LineWidth',2.0, ...
        'DisplayName','Orbital altitude');
    yline(HEIGHT_KM, ':k', 'LineWidth',1.0, 'HandleVisibility','off');
    xlabel('Geodetic latitude (deg)');
    ylabel('Geodetic altitude (km)');
    title(sprintf('Circular orbit: altitude vs latitude   (h = %g km, i = %g°)', ...
        HEIGHT_KM, INC_DEG), 'FontWeight','bold');
    legend('Location','best');
    grid on; xlim([-90 90]);

    fname1 = fullfile(out_dir, sprintf('circular_orbit_altitude_%dkm_i%.0f.png', ...
        HEIGHT_KM, INC_DEG));
    exportgraphics(f1, fname1, 'Resolution',300);
    fprintf('Figure 1 saved -> %s\n', fname1);

    %% ---- Figure 2: altitude + contributing effects on second axis --------
    f2 = figure('Color','w', 'Position',[120 80 820 460]);

    ax_left = axes(f2);
    hold(ax_left, 'on');
    p_alt = plot(ax_left, lat_s, alt_s, '-', 'Color',[0.00 0.45 0.74], ...
        'LineWidth',2.0, 'DisplayName','Orbital altitude');
    yline(ax_left, HEIGHT_KM, ':k', 'LineWidth',1.0, 'HandleVisibility','off');
    ylabel(ax_left, 'Geodetic altitude (km)');
    xlabel(ax_left, 'Geodetic latitude (deg)');
    set(ax_left, 'XColor','k', 'YColor',[0.00 0.45 0.74], 'Box','off');
    grid(ax_left, 'on');
    xlim(ax_left, [-90 90]);

    ax_right = axes(f2, 'Position', ax_left.Position, ...
        'XAxisLocation','top', 'YAxisLocation','right', ...
        'Color','none', 'XTick',[], 'Box','off');
    hold(ax_right, 'on');
    linkaxes([ax_left, ax_right], 'x');

    p_wgs = plot(ax_right, lat_s, wgs84_dev, '-', 'Color',[0.85 0.33 0.10], ...
        'LineWidth',1.6, 'DisplayName','WGS84 surface deviation');
    p_j2  = plot(ax_right, lat_s, j2_dev,    '--','Color',[0.47 0.67 0.19], ...
        'LineWidth',1.6, 'DisplayName','J2 radial perturbation');
    ylabel(ax_right, 'Deviation from mean (km)');
    set(ax_right, 'YColor','k');

    title(ax_left, sprintf('Circular orbit: altitude and contributing effects   (h = %g km, i = %g°)', ...
        HEIGHT_KM, INC_DEG), 'FontWeight','bold');
    legend([p_alt, p_wgs, p_j2], 'Location','best');

    fname2 = fullfile(out_dir, sprintf('circular_orbit_effects_%dkm_i%.0f.png', ...
        HEIGHT_KM, INC_DEG));
    exportgraphics(f2, fname2, 'Resolution',300);
    fprintf('Figure 2 saved -> %s\n', fname2);
end
