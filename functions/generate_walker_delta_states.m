function [r0, v0] = generate_walker_delta_states(orbit_height_m, inclination, total_sats, planes, phasing)
% GENERATE_WALKER_DELTA_STATES  Frozen-orbit ECI initial states for a Walker Delta.
%
%   Standard Walker Delta notation: i:T/P/F
%     T = total_sats, P = planes, F = phasing parameter (integer, 0..P-1)
%
%   RAAN spacing     = 360/P  degrees (uniform)
%   In-plane spacing = 360/(T/P) degrees
%   Cross-plane phase= F * 360/T  degrees between adjacent planes
%
%   Uses the same frozen-orbit Newton solver as generate_walker_star_states
%   so all satellites share the same inertial orbit shape under RK4+J2.
%   The returned (r0, v0) are ready to be passed to PROPAGATE_RK4_J2.
%
%   Inputs
%     orbit_height_m - altitude above WGS-84 equatorial radius (m)
%     inclination    - orbital inclination (deg)
%     total_sats     - total number of satellites (T); must be divisible by planes
%     planes         - number of orbital planes (P)
%     phasing        - phasing factor F (integer, 0..P-1)
%
%   Outputs
%     r0  - Nsat × 3   initial ECI positions  (m)   at epoch
%     v0  - Nsat × 3   initial ECI velocities (m/s) at epoch
%
%   Examples
%     % GPS-like: 24 sats, 6 planes, F=1
%     [r0, v0] = generate_walker_delta_states(20200e3, 55, 24, 6, 1);
%
%     % Starlink shell: 72 sats, 22 planes, F=17
%     [r0, v0] = generate_walker_delta_states(550e3, 53, 66, 6, 2);

    arguments
        orbit_height_m (1,1) double
        inclination    (1,1) double
        total_sats     (1,1) double {mustBeInteger, mustBePositive}
        planes         (1,1) double {mustBeInteger, mustBePositive}
        phasing        (1,1) double {mustBeInteger, mustBeNonnegative}
    end

    if mod(total_sats, planes) ~= 0
        error('total_sats (%d) must be divisible by planes (%d).', total_sats, planes);
    end

    mu   = 3.986004418e14;
    Re_m = 6378.137e3;
    J2   = 1.08262668e-3;

    a              = Re_m + orbit_height_m;
    sats_per_plane = total_sats / planes;

    %% --- Uniform RAAN spacing ---
    delta_raan  = 360 / planes;                   % degrees between planes
    delta_nu    = 360 / sats_per_plane;           % in-plane spacing
    delta_phase = phasing * (360 / total_sats);   % cross-plane phase shift

    Nsat = total_sats;

    %% --- Frozen-orbit reference orbit (RAAN = 0 plane) ---
    [r_ref, v_ref, u_ref] = build_frozen_reference_orbit(a, inclination, mu, J2, Re_m);

    %% --- Per-satellite initial states ---
    r0 = zeros(Nsat, 3);
    v0 = zeros(Nsat, 3);
    k  = 0;
    for p = 1:planes
        raan = (p - 1) * delta_raan;
        cR   = cosd(raan);  sR = sind(raan);
        Rz   = [ cR, -sR, 0;  sR, cR, 0;  0, 0, 1];
        for s = 1:sats_per_plane
            nu = mod((s-1)*delta_nu + (p-1)*delta_phase, 360);
            k  = k + 1;
            [r_pl, v_pl] = sample_reference_at_u(r_ref, v_ref, u_ref, deg2rad(nu));
            r0(k,:) = (Rz * r_pl(:)).';
            v0(k,:) = (Rz * v_pl(:)).';
        end
    end
end

%% =========================================================================
function [r_hist, v_hist, u_hist] = build_frozen_reference_orbit(a, inc_deg, mu, J2, Re)
% Newton-solves (ex, ey) for the frozen orbit then integrates one full period.
    ex = 0;  ey = 0;
    tol_m  = 0.5;
    max_it = 10;
    h_fd   = 5e-4;

    f0 = eval_asymmetry(a, inc_deg, ex, ey, mu, J2, Re);
    for it = 1:max_it
        if max(abs(f0)) < tol_m, break; end
        fx = eval_asymmetry(a, inc_deg, ex + h_fd, ey,        mu, J2, Re);
        fy = eval_asymmetry(a, inc_deg, ex,        ey + h_fd, mu, J2, Re);
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

    e_init  = hypot(ex, ey);
    w_init  = atan2d(ey, ex);
    nu_init = mod(-w_init, 360);
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
    r = r0;  v = v0;
    for i = 2:N
        [r, v] = rk4_step(r, v, dt, mu, J2, Re);
        r_hist(i,:) = r;
        v_hist(i,:) = v;
        u_hist(i)   = arg_of_lat(r, inc_deg);
    end
    u_hist = unwrap(u_hist);
end

%% =========================================================================
function f = eval_asymmetry(a, inc_deg, ex, ey, mu, J2, Re)
    e_init  = hypot(ex, ey);
    w_init  = atan2d(ey, ex);
    nu_init = mod(-w_init, 360);
    [r, v]  = kepler_to_eci(a, e_init, inc_deg, 0, w_init, nu_init, mu);
    r0_mag  = norm(r);

    T  = 2*pi * sqrt(a^3 / mu);
    dt = 1.0;
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
    u_unwrap = unwrap(u_hist);
    r_north  = interp1(u_unwrap, r_mag_hist, pi/2,   'linear', 'extrap');
    r_desc   = interp1(u_unwrap, r_mag_hist, pi,     'linear', 'extrap');
    r_south  = interp1(u_unwrap, r_mag_hist, 3*pi/2, 'linear', 'extrap');
    f = [ r_north - r_south;  r_desc - r0_mag ];
end

%% =========================================================================
function u = arg_of_lat(r, inc_deg)
    ci = cosd(inc_deg);  si = sind(inc_deg);
    x_hat_op = [1, 0, 0];
    y_hat_op = [0, ci, si];
    u = mod(atan2(dot(r, y_hat_op), dot(r, x_hat_op)), 2*pi);
end

%% =========================================================================
function [r_at, v_at] = sample_reference_at_u(r_hist, v_hist, u_hist, u_target)
    if u_target < u_hist(1)
        r_at = r_hist(1,:);  v_at = v_hist(1,:);  return;
    end
    if u_target > u_hist(end)
        r_at = r_hist(end,:);  v_at = v_hist(end,:);  return;
    end
    r_at = [ interp1(u_hist, r_hist(:,1), u_target, 'linear'), ...
             interp1(u_hist, r_hist(:,2), u_target, 'linear'), ...
             interp1(u_hist, r_hist(:,3), u_target, 'linear') ];
    v_at = [ interp1(u_hist, v_hist(:,1), u_target, 'linear'), ...
             interp1(u_hist, v_hist(:,2), u_target, 'linear'), ...
             interp1(u_hist, v_hist(:,3), u_target, 'linear') ];
end

%% =========================================================================
function [r_eci, v_eci] = kepler_to_eci(a, e, inc_deg, raan_deg, w_deg, nu_deg, mu)
    p   = a * (1 - e^2);
    cnu = cosd(nu_deg);  snu = sind(nu_deg);
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

%% =========================================================================
function [r1, v1] = rk4_step(r, v, dt, mu, J2, Re)
    [k1r, k1v] = derivs(r,               v,               mu, J2, Re);
    [k2r, k2v] = derivs(r + 0.5*dt*k1r,  v + 0.5*dt*k1v,  mu, J2, Re);
    [k3r, k3v] = derivs(r + 0.5*dt*k2r,  v + 0.5*dt*k2v,  mu, J2, Re);
    [k4r, k4v] = derivs(r +     dt*k3r,  v +     dt*k3v,  mu, J2, Re);
    r1 = r + (dt/6) * (k1r + 2*k2r + 2*k3r + k4r);
    v1 = v + (dt/6) * (k1v + 2*k2v + 2*k3v + k4v);
end

function [drdt, dvdt] = derivs(r, v, mu, J2, Re)
    rn   = sqrt(sum(r.^2, 2));
    a_tb = -mu .* r ./ rn.^3;
    z2   = r(:,3).^2;
    fac  = -1.5 * J2 * mu * Re^2 ./ rn.^5;
    a_j2 = [ fac .* r(:,1) .* (1 - 5*z2./rn.^2), ...
             fac .* r(:,2) .* (1 - 5*z2./rn.^2), ...
             fac .* r(:,3) .* (3 - 5*z2./rn.^2) ];
    drdt = v;
    dvdt = a_tb + a_j2;
end
