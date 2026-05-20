function [r_ecef_hist, sample_times] = propagate_rk4_j2(r0, v0, start_time, duration_s, dt, sample_step)
% PROPAGATE_RK4_J2  Vectorised fixed-step RK4 propagator with J2 perturbation.
%
%   Takes initial ECI states for any number of satellites and integrates the
%   two-body + J2 equations of motion at a fixed step dt.  At each output
%   sample the ECI position is converted to ECEF via the IAU-1982 GAST
%   approximation, ready for use with ecef2lla.
%
%   Inputs
%     r0          - Nsat × 3   initial ECI positions  (m)
%     v0          - Nsat × 3   initial ECI velocities (m/s)
%     start_time  - datetime   propagation epoch (UTC)
%     duration_s  - scalar     total propagation length (s)
%     dt          - scalar     RK4 step size (s)            [default 1]
%     sample_step - scalar     ECEF output interval (s)     [default 30]
%
%   Outputs
%     r_ecef_hist  - Nsat × 3 × n_samples   ECEF positions (m)
%                   n_samples = floor(duration_s / sample_step) + 1
%     sample_times - 1 × n_samples datetime  UTC epoch of each sample
%
%   Example
%     [r0, v0] = generate_walker_star_states(1000e3, 90, 5, 13, 20);
%     [r_ecef, t] = propagate_rk4_j2(r0, v0, datetime('now','TimeZone','UTC'), 3600);
%     lla = ecef2lla(reshape(permute(r_ecef,[1 3 2]), [], 3));

    arguments
        r0          (:,3) double
        v0          (:,3) double
        start_time  (1,1) datetime
        duration_s  (1,1) double
        dt          (1,1) double = 1.0
        sample_step (1,1) double = 30.0
    end

    mu   = 3.986004418e14;   % m^3/s^2
    Re_m = 6378.137e3;       % WGS-84 equatorial radius (m)
    J2   = 1.08262668e-3;

    Nsat      = size(r0, 1);
    n_samples = floor(duration_s / sample_step) + 1;
    sub_steps = round(sample_step / dt);

    % Pre-compute Julian dates for all output samples (handles month/year rollover)
    sample_times = start_time + seconds((0:n_samples-1) * sample_step);
    utc_mat = [year(sample_times(:)),  month(sample_times(:)),  day(sample_times(:)), ...
               hour(sample_times(:)), minute(sample_times(:)), second(sample_times(:))];
    jd_vec = juliandate(utc_mat);   % n_samples × 1

    r_ecef_hist = zeros(Nsat, 3, n_samples);
    r_ecef_hist(:,:,1) = eci_to_ecef_bulk(r0, jd_vec(1));

    r = r0;  v = v0;
    for sIdx = 2:n_samples
        for k = 1:sub_steps
            [r, v] = rk4_step(r, v, dt, mu, J2, Re_m);
        end
        r_ecef_hist(:,:,sIdx) = eci_to_ecef_bulk(r, jd_vec(sIdx));
    end
end

%% =========================================================================
function [r1, v1] = rk4_step(r, v, dt, mu, J2, Re)
% Vectorised 4th-order Runge-Kutta step.  r, v are Nsat × 3.
    [k1r, k1v] = derivs(r,               v,               mu, J2, Re);
    [k2r, k2v] = derivs(r + 0.5*dt*k1r,  v + 0.5*dt*k1v,  mu, J2, Re);
    [k3r, k3v] = derivs(r + 0.5*dt*k2r,  v + 0.5*dt*k2v,  mu, J2, Re);
    [k4r, k4v] = derivs(r +     dt*k3r,  v +     dt*k3v,  mu, J2, Re);
    r1 = r + (dt/6) * (k1r + 2*k2r + 2*k3r + k4r);
    v1 = v + (dt/6) * (k1v + 2*k2v + 2*k3v + k4v);
end

function [drdt, dvdt] = derivs(r, v, mu, J2, Re)
% Two-body + J2 acceleration, vectorised over Nsat.
    rn   = sqrt(sum(r.^2, 2));        % Nsat × 1
    rn3  = rn.^3;
    rn5  = rn.^5;
    a_tb = -mu .* r ./ rn3;
    z2   = r(:,3).^2;
    fac  = -1.5 * J2 * mu * Re^2 ./ rn5;
    a_j2 = [ fac .* r(:,1) .* (1 - 5*z2./rn.^2), ...
             fac .* r(:,2) .* (1 - 5*z2./rn.^2), ...
             fac .* r(:,3) .* (3 - 5*z2./rn.^2) ];
    drdt = v;
    dvdt = a_tb + a_j2;
end

%% =========================================================================
function r_ecef = eci_to_ecef_bulk(r_eci, jd)
% ECI (J2000) → ECEF via IAU-1982 GAST approximation (~1 arcsec accuracy).
%   r_eci  : Nsat × 3 (m),  jd : scalar Julian date,  r_ecef : Nsat × 3 (m)
    T     = (jd - 2451545.0) / 36525;
    theta = 280.46061837 ...
          + 360.98564736629 * (jd - 2451545.0) ...
          + 0.000387933     * T^2 ...
          - T^3 / 38710000;
    theta  = deg2rad(mod(theta, 360));
    ct = cos(theta);  st = sin(theta);
    R  = [ ct,  st, 0; -st, ct, 0; 0, 0, 1];
    r_ecef = (R * r_eci')';
end
