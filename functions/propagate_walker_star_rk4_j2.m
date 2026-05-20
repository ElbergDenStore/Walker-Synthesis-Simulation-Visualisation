function [all_lat, all_alt_km, t_build, t_prop] = ...
        propagate_walker_star_rk4_j2(orbit_height_m, inclination, planes, sats_per_plane, ...
                                     min_elevation_deg, start_time, duration_s, dt, sample_step)
% PROPAGATE_WALKER_STAR_RK4_J2
%   Builds an asymmetric Walker Star and propagates every satellite with a
%   simple RK4 integrator including the J2 zonal perturbation, working in
%   ECI.  At each sample, converts ECI -> ECEF (Aerospace Toolbox eci2ecef)
%   and then ECEF -> LLA (ecef2lla) so altitude as a function of latitude
%   can be extracted exactly the same way as for the built-in propagators.
%
%   Inputs
%     orbit_height_m     - scalar, metres above equatorial radius
%     inclination        - deg
%     planes             - integer
%     sats_per_plane     - integer
%     min_elevation_deg  - deg (used only for the seam-ratio geometry)
%     start_time         - datetime (UTC)
%     duration_s         - propagation length in seconds
%     dt                 - RK4 step size (s), default 1
%     sample_step        - store output every sample_step seconds, default 30
%
%   Outputs
%     all_lat     - column of geodetic latitudes (deg) across all sats / samples
%     all_alt_km  - column of geodetic altitudes (km), same length
%     t_build     - wall time for initial state setup (s)
%     t_prop      - wall time for RK4 + frame conversion (s)

    arguments
        orbit_height_m     (1,1) double
        inclination        (1,1) double
        planes             (1,1) double {mustBeInteger, mustBePositive}
        sats_per_plane     (1,1) double {mustBeInteger, mustBePositive}
        min_elevation_deg  (1,1) double
        start_time         (1,1) datetime
        duration_s         (1,1) double
        dt                 (1,1) double = 1.0
        sample_step        (1,1) double = 30.0
    end

    %% ------------ Constants ------------------------------------------------
    mu     = 3.986004418e14;     % m^3/s^2
    Re_m   = 6378.137e3;         % WGS-84 equatorial radius (matches eci2ecef)
    J2     = 1.08262668e-3;

    %% ------------ Build initial ECI states ---------------------------------
    tic;
    a = Re_m + orbit_height_m;

    % Seam ratio (identical to the generator)
    Rs_km      = (Re_m + orbit_height_m) / 1e3;
    Re_km      = Re_m / 1e3;
    alpha      = asind((Re_km / Rs_km) * cosd(min_elevation_deg));
    lambda_max = deg2rad(180 - (90 + min_elevation_deg + alpha));
    S          = (2*pi) / sats_per_plane;
    lambda_str = acos(min(1, cos(lambda_max) / cos(S/2)));
    seam_ratio = (2*lambda_str) / (lambda_str + lambda_max);

    co_rotating_spacing = 180 / (planes - 1 + seam_ratio);
    in_plane_spacing    = 360 / sats_per_plane;
    phase_shift         = in_plane_spacing / 2;

    Nsat = planes * sats_per_plane;

    % ----- Per-satellite frozen-orbit osculating initial state -----
    % For each satellite we need (r0, v0) at t=0 such that all satellites
    % share THE SAME inertial orbit (in their respective planes) and that
    % orbit is genuinely frozen under our J2 integrator.  We get this by:
    %   1) Newton-solving a reference (RAAN=0) orbit's initial osculating
    %      (e, omega) so that the orbit is exactly symmetric (N-S and
    %      asc-desc geocentric radii match);
    %   2) Integrating that reference orbit for one full period and storing
    %      (r, v) sampled finely;
    %   3) For each target satellite, looking up the reference state at the
    %      same argument of latitude u_k = nu_k (omega ~ 0 for frozen), then
    %      rotating about ECI z-axis by the satellite's RAAN.
    [r_ref_hist, v_ref_hist, u_ref_hist] = ...
        build_frozen_reference_orbit(a, inclination, mu, J2, Re_m);

    r0 = zeros(Nsat, 3);
    v0 = zeros(Nsat, 3);

    k = 0;
    for p = 1:planes
        raan = (p - 1) * co_rotating_spacing;
        cR   = cosd(raan);  sR = sind(raan);
        Rz   = [ cR, -sR, 0;  sR, cR, 0;  0, 0, 1];   % rotates +x toward +RAAN
        for s = 1:sats_per_plane
            nu = mod((s-1)*in_plane_spacing + (p-1)*phase_shift, 360);
            k  = k + 1;
            [r_pl, v_pl] = sample_reference_at_u(r_ref_hist, v_ref_hist, ...
                                                 u_ref_hist, deg2rad(nu));
            r0(k,:) = (Rz * r_pl(:)).';
            v0(k,:) = (Rz * v_pl(:)).';
        end
    end
    t_build = toc;

    %% ------------ Propagate (vectorised across all sats) ------------------
    tic;
    n_samples = floor(duration_s / sample_step) + 1;
    sub_steps = round(sample_step / dt);

    % Sampled ECI position history: [Nsat x 3 x n_samples]
    r_eci_hist = zeros(Nsat, 3, n_samples);
    r_eci_hist(:,:,1) = r0;

    r = r0;
    v = v0;
    for sIdx = 2:n_samples
        for k_sub = 1:sub_steps
            [r, v] = rk4_step(r, v, dt, mu, J2, Re_m);
        end
        r_eci_hist(:,:,sIdx) = r;
    end

    % Julian dates for each sample (used for GAST rotation)
    sample_times = start_time + seconds((0:n_samples-1) * sample_step);
    utc_mat = [year(sample_times(:)), month(sample_times(:)), day(sample_times(:)), ...
               hour(sample_times(:)), minute(sample_times(:)), second(sample_times(:))];
    jd_vec = juliandate(utc_mat);   % n_samples x 1

    %% ------------ ECI -> ECEF -> LLA --------------------------------------
    all_lat    = zeros(Nsat * n_samples, 1);
    all_alt_km = zeros(Nsat * n_samples, 1);

    for sIdx = 1:n_samples
        r_eci_block = squeeze(r_eci_hist(:,:,sIdx));   % Nsat x 3
        r_ecef      = eci_to_ecef_bulk(r_eci_block, jd_vec(sIdx));  % Nsat x 3
        lla         = ecef2lla(r_ecef);                % Nsat x 3

        rows = (sIdx-1)*Nsat + (1:Nsat);
        all_lat(rows)    = lla(:,1);
        all_alt_km(rows) = lla(:,3) / 1e3;
    end
    t_prop = toc;
end

%% ========================================================================
function [r1, v1] = rk4_step(r, v, dt, mu, J2, Re)
% Vectorised RK4 step.  r and v are Nsat x 3 matrices.
    [k1r, k1v] = derivs(r,              v,              mu, J2, Re);
    [k2r, k2v] = derivs(r + 0.5*dt*k1r, v + 0.5*dt*k1v, mu, J2, Re);
    [k3r, k3v] = derivs(r + 0.5*dt*k2r, v + 0.5*dt*k2v, mu, J2, Re);
    [k4r, k4v] = derivs(r +     dt*k3r, v +     dt*k3v, mu, J2, Re);

    r1 = r + (dt/6) * (k1r + 2*k2r + 2*k3r + k4r);
    v1 = v + (dt/6) * (k1v + 2*k2v + 2*k3v + k4v);
end

function [drdt, dvdt] = derivs(r, v, mu, J2, Re)
% Two-body + J2 acceleration, vectorised.  r, v are Nsat x 3.
    rn   = sqrt(sum(r.^2, 2));       % Nsat x 1
    rn3  = rn.^3;
    rn5  = rn.^5;

    a_tb = -mu .* r ./ rn3;          % two-body, Nsat x 3

    z2   = r(:,3).^2;
    fac  = -1.5 * J2 * mu * Re^2 ./ rn5;
    a_j2 = [ fac .* r(:,1) .* (1 - 5*z2./rn.^2), ...
             fac .* r(:,2) .* (1 - 5*z2./rn.^2), ...
             fac .* r(:,3) .* (3 - 5*z2./rn.^2) ];

    drdt = v;
    dvdt = a_tb + a_j2;
end

%% ========================================================================
%% ========================================================================
function r_ecef = eci_to_ecef_bulk(r_eci_Nx3, jd)
% Vectorised ECI (J2000) -> ECEF via GAST rotation about Z-axis.
% Accurate to ~1 arcsec; no toolbox function needed beyond juliandate.
%   r_eci_Nx3 : Nsat x 3 (metres)
%   jd        : scalar Julian date
%   r_ecef    : Nsat x 3 (metres)
    T     = (jd - 2451545.0) / 36525;          % Julian centuries from J2000
    % IAU 1982 GAST approximation (degrees)
    theta = 280.46061837 ...
          + 360.98564736629 * (jd - 2451545.0) ...
          + 0.000387933     * T^2 ...
          - T^3 / 38710000;
    theta = deg2rad(mod(theta, 360));
    ct = cos(theta);  st = sin(theta);
    % Rz(theta): rotates ECI -> ECEF
    R = [ ct,  st, 0; ...
         -st,  ct, 0; ...
           0,   0, 1];
    r_ecef = (R * r_eci_Nx3')';   % Nsat x 3
end

%% ========================================================================
function [r_eci, v_eci] = kepler_to_eci(a, e, inc_deg, raan_deg, w_deg, nu_deg, mu)
% Classical Keplerian -> ECI (perifocal then 3-1-3 rotation to inertial).
    p   = a * (1 - e^2);
    cnu = cosd(nu_deg);
    snu = sind(nu_deg);
    rpf = [ p*cnu/(1 + e*cnu);  p*snu/(1 + e*cnu);  0 ];
    vpf = sqrt(mu/p) * [ -snu;  e + cnu;  0 ];

    cO = cosd(raan_deg); sO = sind(raan_deg);
    ci = cosd(inc_deg);  si = sind(inc_deg);
    cw = cosd(w_deg);    sw = sind(w_deg);

    R = [ cO*cw - sO*ci*sw,  -cO*sw - sO*ci*cw,   sO*si;
          sO*cw + cO*ci*sw,  -sO*sw + cO*ci*cw,  -cO*si;
                  si*sw,             si*cw,          ci ];

    r_eci = (R * rpf)';
    v_eci = (R * vpf)';
end

%% ========================================================================
function [r_hist, v_hist, u_hist] = build_frozen_reference_orbit(a, inc_deg, mu, J2, Re)
% Builds the reference orbit (RAAN = 0, plane in x-z for i=90 deg).
%
%   Step 1: Newton-solve the osculating eccentricity components (ex, ey) at
%           u = 0 that drive both N-S and asc-desc asymmetry in geocentric
%           radius to zero under our own RK4+J2 integrator.
%   Step 2: integrate one full orbital period with those initial conditions
%           at fine dt (1 s) and return the position/velocity history along
%           with the argument of latitude u at each sample.
%
%   The returned histories can be sampled at any desired u in [0, 2*pi) to
%   produce the t = 0 osculating state of any other satellite in the same
%   plane (just at a different phase along the same inertial orbit).

    % --- 2x2 Newton on (ex, ey) ---
    ex = 0;  ey = 0;
    tol_m  = 0.5;
    max_it = 10;
    h_fd   = 5e-4;

    f0 = eval_asymmetry(a, inc_deg, ex, ey, mu, J2, Re);
    for it = 1:max_it
        if max(abs(f0)) < tol_m, break; end
        fx = eval_asymmetry(a, inc_deg, ex + h_fd, ey, mu, J2, Re);
        fy = eval_asymmetry(a, inc_deg, ex, ey + h_fd, mu, J2, Re);
        Jm = [(fx - f0)/h_fd, (fy - f0)/h_fd];
        if rcond(Jm) < 1e-10, break; end
        d  = -Jm \ f0;
        s  = min(1, 5e-3 / max(abs(d)));
        ex = ex + s*d(1);
        ey = ey + s*d(2);
        f0 = eval_asymmetry(a, inc_deg, ex, ey, mu, J2, Re);
    end
    fprintf('  [frozen-init] Newton converged in %d iter: (e, w) = (%.3e, %+.3f deg)  residual = %.2f m\n', ...
        it, hypot(ex, ey), atan2d(ey, ex), max(abs(f0)));

    % --- Final reference orbit, sampled finely over one period ---
    e_init   = hypot(ex, ey);
    w_init   = atan2d(ey, ex);
    nu_init  = mod(-w_init, 360);     % start exactly at u = w + nu = 0
    [r0, v0] = kepler_to_eci(a, e_init, inc_deg, 0, w_init, nu_init, mu);

    T  = 2*pi * sqrt(a^3 / mu);
    dt = 1.0;
    N  = ceil(T / dt) + 1;
    r_hist = zeros(N, 3);
    v_hist = zeros(N, 3);
    u_hist = zeros(N, 1);
    r_hist(1,:) = r0;
    v_hist(1,:) = v0;
    u_hist(1)   = 0;
    r = r0; v = v0;
    for i = 2:N
        [r, v] = rk4_step(r, v, dt, mu, J2, Re);
        r_hist(i,:) = r;
        v_hist(i,:) = v;
        u_hist(i)   = arg_of_lat(r, inc_deg);
    end
    % Unwrap u so it is monotone in [0, 2*pi) without the atan2 jump.
    u_hist = unwrap(u_hist);
end

%% ========================================================================
function f = eval_asymmetry(a, inc_deg, ex, ey, mu, J2, Re)
% Integrates one orbit with osculating (e, w) = (hypot(ex,ey), atan2d(ey,ex))
% starting at u = 0 (perigee at ascending equator if w=0).  Returns
%   f = [ r(u=pi/2) - r(u=3pi/2);     % N-S radius asymmetry (m)
%         r(u=pi)   - r(u=0)      ]    % asc-desc radius asymmetry (m)
%
% NB: r(u=0) is measured from the IC; r at all other u is recovered by
% interpolation on a u_hist that we UNWRAP to be monotone in [0, 2*pi),
% avoiding the spurious sign change at the natural atan2 wrap.
    e_init  = hypot(ex, ey);
    w_init  = atan2d(ey, ex);
    nu_init = mod(-w_init, 360);
    [r, v]  = kepler_to_eci(a, e_init, inc_deg, 0, w_init, nu_init, mu);
    r0_mag  = norm(r);

    T  = 2*pi * sqrt(a^3 / mu);
    dt = 1.0;                    % match the production-run dt
    N  = ceil(T / dt) + 1;
    r_mag_hist = zeros(N, 1);
    u_hist     = zeros(N, 1);
    r_mag_hist(1) = r0_mag;
    u_hist(1)     = 0;
    for i = 2:N
        [r, v] = rk4_step(r, v, dt, mu, J2, Re);
        r_mag_hist(i) = norm(r);
        u_hist(i)     = arg_of_lat(r, inc_deg);
    end

    % Unwrap u to be monotone increasing from 0 to ~2*pi.  arg_of_lat
    % returns values in [0, 2*pi) so naive plots have a 2*pi jump near the
    % wrap; "unwrap" removes that.
    u_unwrap = unwrap(u_hist);

    r_north = interp1(u_unwrap, r_mag_hist, pi/2,   'linear', 'extrap');
    r_desc  = interp1(u_unwrap, r_mag_hist, pi,     'linear', 'extrap');
    r_south = interp1(u_unwrap, r_mag_hist, 3*pi/2, 'linear', 'extrap');

    f = [ r_north - r_south;
          r_desc  - r0_mag ];
end

%% ========================================================================
function u = arg_of_lat(r, inc_deg)
% Argument of latitude for an orbit with RAAN = 0.  For i = 90 deg the
% orbit plane is the x-z plane, ascending node along +x, north (u = pi/2)
% along +z.  Generalised for any i by projecting onto the orbital frame.
    % Orbital normal for RAAN=0: n = (sin i)*y_hat ... actually
    %   n = Rz(RAAN) * Rx(i) * z_hat = (0, -sin i, cos i)   (for RAAN=0)
    % "Ascending node direction" = +x.  In-plane "north" direction =
    %   Rx(i) * y_hat = (0, cos i, sin i).
    ci = cosd(inc_deg);  si = sind(inc_deg);
    x_hat_op = [1, 0, 0];                % ascending node direction
    y_hat_op = [0, ci, si];              % 90 deg along the orbit from node
    a1 = dot(r, x_hat_op);
    a2 = dot(r, y_hat_op);
    u  = mod(atan2(a2, a1), 2*pi);
end

%% ========================================================================
function [r_at, v_at] = sample_reference_at_u(r_hist, v_hist, u_hist, u_target)
% Returns reference-orbit (r, v) at the requested argument of latitude.
% u_hist is assumed already unwrapped (monotone in [0, 2*pi)).  Per-component
% linear interpolation; for u_target = 0 we fall back to the first sample
% (the IC), for u_target = 2*pi we extrapolate to the end of the orbit.
    if u_target < u_hist(1)
        r_at = r_hist(1,:);
        v_at = v_hist(1,:);
        return;
    end
    if u_target > u_hist(end)
        r_at = r_hist(end,:);
        v_at = v_hist(end,:);
        return;
    end
    r_at = [ interp1(u_hist, r_hist(:,1), u_target, 'linear'), ...
             interp1(u_hist, r_hist(:,2), u_target, 'linear'), ...
             interp1(u_hist, r_hist(:,3), u_target, 'linear') ];
    v_at = [ interp1(u_hist, v_hist(:,1), u_target, 'linear'), ...
             interp1(u_hist, v_hist(:,2), u_target, 'linear'), ...
             interp1(u_hist, v_hist(:,3), u_target, 'linear') ];
end
