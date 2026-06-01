function results = validate_satellite_initialization(HEIGHT_KM, INC_DEG)
% VALIDATE_SATELLITE_INITIALIZATION
%   Shows how non-circular the propagated orbit really is, by comparing the
%   altitude on the ascending vs descending equator crossing, and on the
%   north vs south pole crossing.
%
%   Propagators compared:
%       sgp4                  - treats inputs as Brouwer-Lyddane mean
%                               elements; e=0,w=0 -> residual eccentricity
%                               vector along the POLAR axis -> N-S asymmetry.
%       numerical             - integrates the input osculating state.
%                               e=0,w=0 puts radial velocity = 0 at the
%                               equator, which is NOT what a J2-mean circle
%                               looks like there -> residual eccentricity
%                               along the EQUATORIAL axis -> asc-desc
%                               asymmetry.
%       two-body-keplerian    - no J2, r = const exactly. Reference only.
%       numerical-frozen      - same numerical integrator, but the initial
%                               osculating eccentricity vector (ex, ey) is
%                               solved by 2x2 Newton so that both N-S and
%                               asc-desc asymmetries are driven to zero.
%
%   Usage:
%       validate_satellite_initialization()           % 1000 km, 90 deg
%       validate_satellite_initialization(600, 87.9)

    if nargin < 1 || isempty(HEIGHT_KM), HEIGHT_KM = 1000; end
    if nargin < 2 || isempty(INC_DEG),   INC_DEG   = 90;   end

    Re_eq_m = 6378.14e3;
    mu      = 3.986004418e14;
    a_m     = Re_eq_m + HEIGHT_KM*1e3;
    T_s     = 2*pi*sqrt(a_m^3 / mu);

    sc = satelliteScenario;
    sc.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone','UTC');
    sc.StopTime   = sc.StartTime + seconds(T_s*2 + 5);   % two orbits
    sc.SampleTime = 5;

    fprintf('\nRequested: a = R_eq + h = %.3f km   (h = %g km, i = %g deg)\n', ...
        a_m/1e3, HEIGHT_KM, INC_DEG);
    fprintf('Period    : %.2f min\n\n', T_s/60);

    %% Baseline: three propagators with e=0, w=0
    propagators = {'sgp4', 'numerical', 'two-body-keplerian'};
    field_names = {'sgp4', 'numerical', 'keplerian'};

    results = struct();
    for p = 1:numel(propagators)
        s = run_one(sc, a_m, 0, INC_DEG, 0, 0, propagators{p});
        s.label = propagators{p};
        results.(field_names{p}) = s;
    end

    %% Newton-solve the numerical propagator for a J2-mean circular orbit
    [ex_star, ey_star, n_iter] = solve_frozen(sc, a_m, INC_DEG);
    e_star = hypot(ex_star, ey_star);
    w_star = atan2d(ey_star, ex_star);
    s_fr = run_one(sc, a_m, e_star, INC_DEG, 0, w_star, 'numerical');
    s_fr.label = 'numerical-frozen';
    results.numerical_frozen = s_fr;

    %% Print headline table
    order = {'sgp4','numerical','keplerian','numerical_frozen'};
    fprintf('  %-20s   r_pp [m]   alt_pp [km]   N-S [m]   asc-desc [m]   e_osc\n', 'propagator');
    fprintf('  ------------------------------------------------------------------------------\n');
    for k = 1:numel(order)
        s = results.(order{k});
        fprintf('  %-20s  %8.1f    %8.3f    %+8.1f    %+8.1f      %.2e\n', ...
            s.label, s.r_pp_m, s.alt_pp_km, s.ns_m, s.ad_m, s.e_osc);
    end
    fprintf('\n  r_pp     = peak-peak geocentric radius (0 = perfectly circular)\n');
    fprintf('  alt_pp   = peak-peak geodetic altitude (~21.4 km is the WGS84 floor)\n');
    fprintf('  N-S      = alt(north pole) - alt(south pole)\n');
    fprintf('  asc-desc = alt(equator, ascending) - alt(equator, descending)\n');
    fprintf('  e_osc    = (r_max - r_min) / (r_max + r_min)\n\n');

    fprintf('Frozen-orbit solver: %d Newton iters\n', n_iter);
    fprintf('  e_init = %.3e   w_init = %+7.2f deg   (nu_init = -w)\n\n', e_star, w_star);

    %% Plot - SGP4 vs numerical vs frozen-numerical
    f = figure('Color','w', 'Position',[80 100 980 560]);
    hold on;
    plot_traces(results.sgp4,             [0.85 0.33 0.10], 'SGP4');
    plot_traces(results.numerical,        [0.47 0.67 0.19], 'numerical (e=0)');
    plot_traces(results.numerical_frozen, [0.00 0.45 0.74], 'numerical frozen');
    yline(HEIGHT_KM, ':k', 'LineWidth',1.0, 'HandleVisibility','off');
    xlabel('Geodetic latitude (deg)');
    ylabel('Geodetic altitude (km)');
    title(sprintf('Altitude vs latitude   h_{req}=%g km, i=%g deg', HEIGHT_KM, INC_DEG), ...
        'FontWeight','bold');
    legend('Location','best');
    grid on; xlim([-90 90]);

    script_dir = fileparts(mfilename('fullpath'));
    out_dir    = fullfile(script_dir, 'figures');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    fname = fullfile(out_dir, sprintf('eccentricity_check_%dkm_i%.0f.png', HEIGHT_KM, INC_DEG));
    exportgraphics(f, fname, 'Resolution',200);
    fprintf('Plot saved -> %s\n', fname);
end

%% --------------------------------------------------------------------
function s = run_one(sc, a, e, inc, raan, w_deg, prop)
% Build one sat, propagate, return analysis struct + raw traces.
%   nu0 is chosen so the satellite sits at u = w + nu = 0  (ascending eq).
    nu_deg = mod(-w_deg, 360);
    sat = satellite(sc, a, e, inc, raan, w_deg, nu_deg, ...
        'Name', sprintf('%s_e%.1e_w%+05.1f', prop, e, w_deg), ...
        'OrbitPropagator', prop);
    P   = squeeze(states(sat, 'CoordinateFrame','ECEF'))';
    lla = ecef2lla(P);
    s   = analyse(lla(:,1), lla(:,3)/1e3, sqrt(sum(P.^2,2))/1e3);
    s.lat = lla(:,1); s.alt = lla(:,3)/1e3; s.r = sqrt(sum(P.^2,2))/1e3;
end

%% --------------------------------------------------------------------
function [ex, ey, it] = solve_frozen(sc, a, inc)
% 2x2 Newton on (ex, ey) -> ([ns_m; ad_m] = 0).
%   Forward differences for the Jacobian, start from (0,0).
    ex = 0; ey = 0;
    tol_m  = 0.5;             % half a metre on both metrics
    max_it = 6;
    h      = 5e-4;            % FD step in eccentricity components
    for it = 1:max_it
        f0 = eval_asym(sc, a, inc, ex, ey);
        if max(abs(f0)) < tol_m, return; end
        fx = eval_asym(sc, a, inc, ex + h, ey);
        fy = eval_asym(sc, a, inc, ex,     ey + h);
        J  = [(fx - f0)/h, (fy - f0)/h];
        if rcond(J) < 1e-10
            return;
        end
        d  = -J \ f0;
        % Trust-region cap: don't let |d| exceed 5e-3 in eccentricity-vector L-inf
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

%% --------------------------------------------------------------------
function plot_traces(s, col, name)
    asc = [false; diff(s.lat) >= 0];

    % Sort by latitude before plotting so MATLAB connects points in spatial
    % order, not time order.  Without this, the two polar regions (Q1 top
    % and Q4 top) are connected by a long horizontal line that visually
    % dominates the figure — especially when the orbit is symmetric and
    % both poles have the same altitude.
    [la, ia] = sort(s.lat( asc)); aa = s.alt( asc); aa = aa(ia);
    [ld, id] = sort(s.lat(~asc)); ad = s.alt(~asc); ad = ad(id);

    plot(la, aa, '-',  'Color',col, 'LineWidth',1.7, ...
        'DisplayName', sprintf('%s ascending',  name));
    plot(ld, ad, '--', 'Color',col, 'LineWidth',1.7, ...
        'DisplayName', sprintf('%s descending', name));
end

%% --------------------------------------------------------------------
function s = analyse(lat, alt, r)
    s.r_pp_m    = (max(r)    - min(r)) * 1e3;
    s.alt_pp_km =  max(alt)  - min(alt);
    s.e_osc     = (max(r) - min(r)) / (max(r) + min(r));

    n_mask = lat >  85;
    s_mask = lat < -85;
    if any(n_mask) && any(s_mask)
        s.ns_m = (mean(alt(n_mask)) - mean(alt(s_mask))) * 1e3;
    else
        s.ns_m = NaN;
    end

    dlat    = [0; diff(lat)];
    eq_mask = abs(lat) < 5;
    asc = eq_mask & dlat > 0;
    des = eq_mask & dlat < 0;
    if any(asc) && any(des)
        s.ad_m = (mean(alt(asc)) - mean(alt(des))) * 1e3;
    else
        s.ad_m = NaN;
    end
end
