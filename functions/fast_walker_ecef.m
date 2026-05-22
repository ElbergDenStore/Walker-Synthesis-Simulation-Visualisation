function sat_pos_ecef = fast_walker_ecef(Orbit_height, Inc_deg, P, S, F, time_steps_sec, StartTime)
% FAST_WALKER_ECEF  Pure-math ECEF propagator for a Walker Delta constellation.
%
%   sat_pos_ecef = fast_walker_ecef(Orbit_height, Inc_deg, P, S, F,
%                                   time_steps_sec, StartTime)
%
%   Produces circular two-body Walker Delta orbits in ECEF without requiring
%   the Satellite Communications Toolbox.
%
%   Inputs
%     Orbit_height    Orbital altitude above Earth surface (m)
%     Inc_deg         Inclination (deg)
%     P               Number of orbital planes
%     S               Number of satellites per plane
%     F               Phasing factor  (0 <= F < P)
%     time_steps_sec  Row vector of simulation times from epoch (s)
%     StartTime       datetime scalar (UTC) used to compute GMST
%
%   Output
%     sat_pos_ecef    [3 x (P*S) x nT] ECEF position matrix (m)

    T = P * S;
    a  = 6378.137e3 + Orbit_height; % Semi-major axis (m)
    mu = 3.986004418e14;            % Earth gravitational constant (m^3/s^2)
    we = 7.2921150e-5;              % Earth rotation rate (rad/s)
    inc = deg2rad(Inc_deg);

    n  = sqrt(mu / a^3);            % Mean motion (rad/s)
    nT = length(time_steps_sec);

    % GMST at StartTime
    JD       = juliandate(StartTime);
    D        = JD - 2451545.0;      % Days since J2000.0
    GMST_deg = mod(280.46061837 + 360.98564736629 * D, 360);
    theta_g0 = deg2rad(GMST_deg);

    sat_pos_ecef = zeros(3, T, nT);

    sat_idx = 1;
    for p = 0:(P-1)
        RAAN = p * (2*pi / P);

        for s = 0:(S-1)
            % Initial Mean Anomaly (= true anomaly for circular orbit)
            M0 = s * (2*pi / S) + p * F * (2*pi / T);

            theta = M0 + n * time_steps_sec;

            x_orb = a * cos(theta);
            y_orb = a * sin(theta);

            X_eci = x_orb * cos(RAAN) - y_orb * cos(inc) * sin(RAAN);
            Y_eci = x_orb * sin(RAAN) + y_orb * cos(inc) * cos(RAAN);
            Z_eci = y_orb * sin(inc);

            theta_g = theta_g0 + we * time_steps_sec;
            sat_pos_ecef(1, sat_idx, :) =  X_eci .* cos(theta_g) + Y_eci .* sin(theta_g);
            sat_pos_ecef(2, sat_idx, :) = -X_eci .* sin(theta_g) + Y_eci .* cos(theta_g);
            sat_pos_ecef(3, sat_idx, :) =  Z_eci;

            sat_idx = sat_idx + 1;
        end
    end
end
