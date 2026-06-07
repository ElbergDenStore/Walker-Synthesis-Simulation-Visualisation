% PLOT_CIRCULAR_ORBIT_ALTITUDE
%   Produces TWO publication-quality figures:
%
%   Figure 1 – Orbital altitude vs geodetic latitude (single pass).
%
%   Figure 2 – Same altitude trace on the primary y-axis, plus a secondary
%   y-axis showing the WGS84 surface deviation: how much the reference
%   ellipsoid radius varies with latitude (Earth's equatorial bulge).
%
%   Propagated with MATLAB's two-body-keplerian propagator (no J2, no drag).
%   With e=0 the orbit is perfectly circular in ECI; all altitude variation
%   seen vs latitude is purely due to the WGS84 oblate ellipsoid.
%
%   Edit HEIGHT_KM and INC_DEG below.

    HEIGHT_KM = 1000;   % orbital altitude (km)
    INC_DEG   = 90;     % inclination (deg)

    %% Constants ------------------------------------------------------------
    Re_m   = 6378.137e3;   % WGS-84 equatorial radius (m)
    b_WGS  = 6356752.3142; % WGS-84 polar radius (m)
    mu     = 3.986004418e14;

    a_m  = Re_m + HEIGHT_KM*1e3;
    T_s  = 2*pi * sqrt(a_m^3 / mu);

    e     = 0;   % perfectly circular
    w_deg = 0;
    nu_deg = 0;

    fprintf('\nCircular-orbit altitude vs latitude\n');
    fprintf('  h = %g km   i = %g deg\n', HEIGHT_KM, INC_DEG);
    fprintf('  e = 0 (two-body-keplerian)   Period: %.2f min\n\n', T_s/60);

    %% Propagate ------------------------------------------------------------
    sc = satelliteScenario;
    sc.StartTime  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
    sc.StopTime   = sc.StartTime + seconds(T_s + 5);   % one orbit
    sc.SampleTime = 5;

    sat = satellite(sc, a_m, e, INC_DEG, 0, w_deg, nu_deg, ...
        'Name', 'circular', 'OrbitPropagator', 'two-body-keplerian');

    P   = squeeze(states(sat, 'CoordinateFrame','ECEF'))';  % Nx3, metres
    lla = ecef2lla(P);
    lat = lla(:,1);          % geodetic latitude (deg)
    alt = lla(:,3) / 1e3;   % geodetic altitude (km)

    % One clean pass: use ascending half, sorted by latitude
    asc_mask    = [false; diff(lat) >= 0];
    [lat_s, ix] = sort(lat(asc_mask));
    tmp         = alt(asc_mask);  alt_s = tmp(ix);

    %% Derived quantities for Figure 2 -------------------------------------
    % WGS84 geocentric surface radius at each latitude (analytical).
    % With a perfectly circular orbit (r=a=const in ECI), the geodetic
    % altitude variation is entirely caused by the oblate WGS84 ellipsoid.
    lat_rad    = deg2rad(lat_s);
    R_wgs84_km = sqrt( ((Re_m^2 * cos(lat_rad)).^2 + (b_WGS^2 * sin(lat_rad)).^2) ./ ...
                       ((Re_m   * cos(lat_rad)).^2 + (b_WGS   * sin(lat_rad)).^2) ) / 1e3;
    wgs84_dev  = R_wgs84_km - mean(R_wgs84_km);   % deviation from its own mean (km)

    %% ---- Figure 1: simple altitude trace --------------------------------
    script_dir = fileparts(mfilename('fullpath'));
    out_dir    = fullfile(script_dir, 'figures');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    f1 = figure('Color','w', 'Position',[80 120 500 250]);
    hold on;
    plot(lat_s, alt_s, '-', 'Color',[0.00 0.45 0.74], 'LineWidth',2.0);
    yline(HEIGHT_KM, ':k', 'LineWidth',1.0, 'HandleVisibility','off');
    xlabel('Geodetic latitude (deg)');
    ylabel('Geodetic altitude (km)');
    title(sprintf('Constant Equatorial Orbital Altitude | h = %g km, i = %g°', ...
        HEIGHT_KM, INC_DEG), 'FontWeight','bold');
    % legend('Location','best');
    grid on; xlim([-90 90]);

    fname1 = fullfile(out_dir, sprintf('circular_orbit_altitude_%dkm_i%.0f.png', ...
        HEIGHT_KM, INC_DEG));
    exportgraphics(f1, fname1, 'Resolution',600);
    fprintf('Figure 1 saved -> %s\n', fname1);

    %% ---- Figure 2: altitude + contributing effects on second axis --------
    f2 = figure('Color','w', 'Position',[80 120 500 250]);

    ax_left = axes(f2);
    hold(ax_left, 'on');
    p_alt = plot(ax_left, lat_s, alt_s, '-', 'Color',[0.00 0.45 0.74], ...
        'LineWidth',2.0,'DisplayName','Circular orbit');
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
    ylabel(ax_right, 'WGS84 deviation from mean (km)');
    set(ax_right, 'YColor',[0.85 0.33 0.10]);

    title(ax_left, sprintf('Circular orbit: altitude and WGS84 effect   (h = %g km, i = %g°)', ...
        HEIGHT_KM, INC_DEG), 'FontWeight','bold');
    legend([p_alt, p_wgs], 'Location','best');

    fname2 = fullfile(out_dir, sprintf('circular_orbit_effects_%dkm_i%.0f.png', ...
        HEIGHT_KM, INC_DEG));
    exportgraphics(f2, fname2, 'Resolution',300);
    fprintf('Figure 2 saved -> %s\n', fname2);
