function sat_pos_ecef = fast_walker_star_ecef(Orbit_height, Inc_deg, P, S, min_elevation_deg, min_lat_deg, time_steps_sec, StartTime)
% FAST_WALKER_STAR_ECEF  Pure-math ECEF propagator for an asymmetrical Walker Star.
%
%   sat_pos_ecef = fast_walker_star_ecef(Orbit_height, Inc_deg, P, S,
%                       min_elevation_deg, min_lat_deg, time_steps_sec, StartTime)
%
%   Produces the same geometry as asymmetrical_walker_star_generation with
%   'two-body-keplerian', without requiring the Satellite Communications
%   Toolbox.
%
%   Inputs
%     Orbit_height      Orbital altitude above Earth surface (m)
%     Inc_deg           Inclination (deg)
%     P                 Number of orbital planes
%     S                 Number of satellites per plane
%     min_elevation_deg Minimum service elevation angle (deg); used for seam ratio
%     min_lat_deg       Coverage reference latitude (deg); 0 = equatorial
%     time_steps_sec    Row vector of simulation times from epoch (s)
%     StartTime         datetime scalar (UTC) used to compute GMST
%
%   Output
%     sat_pos_ecef      [3 x (P*S) x nT] ECEF position matrix (m)

    T = P * S;
    r_earth = 6378.14e3;
    a   = r_earth + Orbit_height;   % Semi-major axis (m)
    mu  = 3.986004418e14;           % Earth gravitational constant (m^3/s^2)
    we  = 7.2921150e-5;             % Earth rotation rate (rad/s)
    inc = deg2rad(Inc_deg);

    n  = sqrt(mu / a^3);            % Mean motion (rad/s)
    nT = length(time_steps_sec);

    % GMST at StartTime
    JD       = juliandate(StartTime);
    D        = JD - 2451545.0;      % Days since J2000.0
    GMST_deg = mod(280.46061837 + 360.98564736629 * D, 360);
    theta_g0 = deg2rad(GMST_deg);

    % --- Seam-ratio geometry (identical to asymmetrical_walker_star_generation) ---
    orbit_height_km = Orbit_height / 1000;
    Re_eq_km = 6378.14;
    Rs_km    = Re_eq_km + orbit_height_km;
    a_wgs = 6378.137;  b_wgs = 6356.7523142;
    lat_r = deg2rad(min_lat_deg);
    Re_km = sqrt((a_wgs^4*cos(lat_r)^2 + b_wgs^4*sin(lat_r)^2) / ...
                 (a_wgs^2*cos(lat_r)^2 + b_wgs^2*sin(lat_r)^2));

    alpha       = asind((Re_km / Rs_km) * cosd(min_elevation_deg));
    lambda_max  = deg2rad(180 - (90 + min_elevation_deg + alpha));
    S_ang       = (2*pi) / S;
    lambda_str  = acos(min(1, cos(lambda_max) / cos(S_ang / 2)));
    seam_ratio  = (2 * lambda_str) / (lambda_str + lambda_max);

    % RAAN spacing and in-plane phase stagger
    co_rot_spacing   = 180 / (P - 1 + seam_ratio);   % degrees between planes
    in_plane_spacing = 360 / S;                        % degrees between sats in plane
    phase_shift      = in_plane_spacing / 2;           % brick-wall stagger (deg)

    sat_pos_ecef = zeros(3, T, nT);

    sat_idx = 1;
    for p = 0:(P-1)
        RAAN = deg2rad(p * co_rot_spacing);

        for s = 0:(S-1)
            nu0   = deg2rad(mod(s * in_plane_spacing + p * phase_shift, 360));
            theta  = nu0 + n * time_steps_sec;

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
