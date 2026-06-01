function sat_pos_ecef = fast_walker_ecef(Orbit_height, Inc_deg, P, S, F, ...
        time_steps_sec, StartTime, WalkerStar, min_elevation_deg, min_lat_deg)
% FAST_WALKER_ECEF  Pure-math ECEF propagator for Walker Delta or Walker Star.
%
%   sat_pos_ecef = fast_walker_ecef(Orbit_height, Inc_deg, P, S, F,
%                                   time_steps_sec, StartTime)
%   sat_pos_ecef = fast_walker_ecef(..., StartTime, WalkerStar,
%                                   min_elevation_deg, min_lat_deg)
%
%   Produces circular two-body orbits in ECEF without requiring the Satellite
%   Communications Toolbox.  A Walker Delta (default) uses uniform RAAN spacing
%   and an integer phasing factor F.  A Walker Star (WalkerStar = true) uses the
%   asymmetric seam-ratio RAAN spacing and brick-wall in-plane stagger, matching
%   generate_walker_star_scenario with 'two-body-keplerian'.
%
%   Inputs
%     Orbit_height      Orbital altitude above Earth surface (m)
%     Inc_deg           Inclination (deg)
%     P                 Number of orbital planes
%     S                 Number of satellites per plane
%     F                 Phasing factor (0 <= F < P).  Delta only; ignored for star.
%     time_steps_sec    Row vector of simulation times from epoch (s)
%     StartTime         datetime scalar (UTC) used to compute GMST
%     WalkerStar        (optional) logical; true for Walker Star. Default false.
%     min_elevation_deg (optional, star only) Minimum service elevation (deg)
%     min_lat_deg       (optional, star only) Coverage reference latitude (deg)
%
%   Output
%     sat_pos_ecef      [3 x (P*S) x nT] ECEF position matrix (m)

    if nargin < 8 || isempty(WalkerStar);        WalkerStar = false; end
    if nargin < 9 || isempty(min_elevation_deg); min_elevation_deg = 0; end
    if nargin < 10 || isempty(min_lat_deg);      min_lat_deg = 0; end

    T   = P * S;
    mu  = 3.986004418e14;           % Earth gravitational constant (m^3/s^2)
    we  = 7.2921150e-5;             % Earth rotation rate (rad/s)
    inc = deg2rad(Inc_deg);

    r_earth = 6378.137e3;           % Earth equatorial radius (m), WGS84
    a  = r_earth + Orbit_height;    % Semi-major axis (m)
    n  = sqrt(mu / a^3);            % Mean motion (rad/s)
    nT = length(time_steps_sec);

    % GMST at StartTime
    JD       = juliandate(StartTime);
    D        = JD - 2451545.0;      % Days since J2000.0
    GMST_deg = mod(280.46061837 + 360.98564736629 * D, 360);
    theta_g0 = deg2rad(GMST_deg);

    % --- Per-plane RAAN and in-plane phase stagger ---
    if WalkerStar
        % Seam-ratio geometry (identical to generate_walker_star_scenario)
        orbit_height_km = Orbit_height / 1000;
        Re_eq_km = 6378.137;
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

        co_rot_spacing   = 180 / (P - 1 + seam_ratio);   % deg between planes
        in_plane_spacing = 360 / S;                       % deg between sats in plane
        phase_shift      = in_plane_spacing / 2;          % brick-wall stagger (deg)
    end

    sat_pos_ecef = zeros(3, T, nT);

    sat_idx = 1;
    for p = 0:(P-1)
        if WalkerStar
            RAAN = deg2rad(p * co_rot_spacing);
        else
            RAAN = p * (2*pi / P);
        end

        for s = 0:(S-1)
            if WalkerStar
                nu0 = deg2rad(mod(s * in_plane_spacing + p * phase_shift, 360));
            else
                % Initial Mean Anomaly (= true anomaly for circular orbit)
                nu0 = s * (2*pi / S) + p * F * (2*pi / T);
            end

            theta = nu0 + n * time_steps_sec;

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
