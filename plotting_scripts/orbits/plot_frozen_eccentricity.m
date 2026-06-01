function plot_frozen_eccentricity()
% SURVEY_FROZEN_ECCENTRICITY
%   Runs the Newton frozen-orbit solver at altitudes from 500 to 1200 km
%   (i = 90 deg) and plots the required initialization eccentricity and
%   argument of perigee as functions of altitude.
%
%   Also overlays the first-order analytical approximation
%       e_frozen ≈ (J2 / 2) * (Re / a)^2
%   so you can judge whether a fixed formula is good enough for your
%   simulator (avoiding the Newton solve at run-time).
%
%   Usage:
%       survey_frozen_eccentricity()

    heights_km = 500 : 50 : 1200;   % 15 points
    nH         = numel(heights_km);

    Re_eq_m = 6378.137e3;
    mu      = 3.986004418e14;
    J2      = 1.0826e-3;

    e_sol = nan(1, nH);
    w_sol = nan(1, nH);

    fprintf('\n%-8s  %-12s  %-10s  iters\n', 'h [km]', 'e_frozen', 'w [deg]');
    fprintf('%s\n', repmat('-', 1, 42));

    for k = 1:nH
        h_km = heights_km(k);
        a_m  = Re_eq_m + h_km * 1e3;
        T_s  = 2*pi * sqrt(a_m^3 / mu);

        sc = satelliteScenario;
        sc.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone','UTC');
        sc.StopTime   = sc.StartTime + seconds(T_s*2 + 5);
        sc.SampleTime = 5;

        [ex, ey, it] = solve_frozen(sc, a_m, 90);
        e_sol(k) = hypot(ex, ey);
        w_sol(k) = atan2d(ey, ex);

        fprintf('%-8g  %-12.4e  %+8.2f       %d\n', h_km, e_sol(k), w_sol(k), it);
    end

    %% Analytical approximation: e ≈ (J2/2)*(Re/a)^2
    a_km  = Re_eq_m/1e3 + heights_km;
    e_ana = (J2/2) * (Re_eq_m/1e3 ./ a_km).^2;

    fprintf('\nRMS error of analytical fit: %.2e\n', rms(e_sol - e_ana));

    %% Plot
    f = figure('Color','w', 'Position',[80 80 820 580]);

    ax1 = subplot(2,1,1); hold on;
    plot(heights_km, e_sol*1e4, 'bo-', 'LineWidth',1.8, ...
        'MarkerFaceColor','b', 'MarkerSize',6, 'DisplayName','Newton solver');
    plot(heights_km, e_ana*1e4, 'r--', 'LineWidth',1.5, ...
        'DisplayName','(J_2/2)(R_e/a)^2  analytical');
    ylabel('e_{frozen}   (\times10^{-4})');
    title('Frozen-orbit initialization eccentricity vs altitude  (i = 90°)', ...
        'FontWeight','bold');
    leg = legend('Location','northeast');
    leg.Interpreter = 'tex';
    grid on;

    ax2 = subplot(2,1,2); hold on;
    plot(heights_km, w_sol, 'bs-', 'LineWidth',1.8, ...
        'MarkerFaceColor','b', 'MarkerSize',6);
    yline(0, ':k', 'LineWidth',1.0);
    xlabel('Orbital altitude (km)');
    ylabel('\omega_0  (deg)');
    title('Argument of perigee at t_0', 'FontWeight','bold');
    grid on;
    % Reasonable y-range; if w drifts far from 0, expand automatically
    yl = max(abs(w_sol)) * 1.5 + 1;
    ylim([-yl yl]);

    linkaxes([ax1 ax2], 'x');

    script_dir = fileparts(mfilename('fullpath'));
    out_dir    = fullfile(script_dir, 'figures');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    fname = fullfile(out_dir, 'frozen_eccentricity_survey.png');
    exportgraphics(f, fname, 'Resolution',200);
    fprintf('Plot saved -> %s\n', fname);
end

%% --------------------------------------------------------------------
function [ex, ey, it] = solve_frozen(sc, a, inc)
    ex = 0; ey = 0;
    tol_m  = 0.5;
    max_it = 6;
    h      = 5e-4;
    for it = 1:max_it
        f0 = eval_asym(sc, a, inc, ex, ey);
        if max(abs(f0)) < tol_m, return; end
        fx = eval_asym(sc, a, inc, ex + h, ey);
        fy = eval_asym(sc, a, inc, ex,     ey + h);
        J  = [(fx - f0)/h, (fy - f0)/h];
        if rcond(J) < 1e-10, return; end
        d  = -J \ f0;
        scale = min(1, 5e-3 / max(abs(d)));
        ex = ex + scale*d(1);
        ey = ey + scale*d(2);
    end
end

function f = eval_asym(sc, a, inc, ex, ey)
    e = hypot(ex, ey);
    w = atan2d(ey, ex);
    s = run_one(sc, a, e, inc, 0, w, 'numerical');
    f = [s.ns_m; s.ad_m];
end

function s = run_one(sc, a, e, inc, raan, w_deg, prop)
    nu_deg = mod(-w_deg, 360);
    sat = satellite(sc, a, e, inc, raan, w_deg, nu_deg, ...
        'Name', sprintf('%s_e%.1e_w%+05.1f', prop, e, w_deg), ...
        'OrbitPropagator', prop);
    P   = squeeze(states(sat, 'CoordinateFrame','ECEF'))';
    lla = ecef2lla(P);
    s.lat   = lla(:,1);
    s.alt   = lla(:,3)/1e3;
    s.r     = sqrt(sum(P.^2,2))/1e3;
    n_mask  = s.lat >  85;
    s_mask  = s.lat < -85;
    s.ns_m  = (mean(s.alt(n_mask)) - mean(s.alt(s_mask))) * 1e3;
    dlat    = [0; diff(s.lat)];
    eq      = abs(s.lat) < 5;
    s.ad_m  = (mean(s.alt(eq & dlat > 0)) - mean(s.alt(eq & dlat < 0))) * 1e3;
end
