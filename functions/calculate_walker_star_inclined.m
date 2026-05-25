function [optimal_planes, optimal_sats_per_plane, total_sats] = ...
    calculate_walker_star_inclined(orbit_height_km, min_latitude_deg, ...
                                   min_elevation_deg, inclination_deg)
% CALCULATE_WALKER_STAR_INCLINED  Minimum-satellite Walker Star for inclination < 90°.
%
%   Uses the spherical street-of-coverage model (Rider/Ballard/Walker) extended
%   for inclined orbits.  For i = 90° the result equals calculate_walker_star().
%
%   Physics:
%     At latitude φ a satellite track drifts in longitude (from Napier's Circle):
%         Δλ(φ,i) = arcsin( tan(φ) / tan(i) )
%     Counter-rotating planes at the seam share an ascending-node RAAN.  Their
%     tracks bow outward by ±Δλ, widening the coverage gap at the seam.
%
%     Seam RAAN budget (spherical geometry at critical latitude φ):
%         R_seam = arccos( (cos(D_ctr) − sin²φ) / cos²φ )
%     Reduction due to inclination (track bowing):
%         R_seam_max = R_seam − 2·Δλ(φ,i)          [must be > 0]
%     Co-rotating RAAN budget per interval:
%         R_co_max = arccos( (cos(D_same) − sin²φ) / cos²φ )
%     Number of planes:
%         P = ceil( (π − R_seam_max) / R_co_max ) + 1
%
%   Returns [Inf, 0, Inf] if inclination_deg <= min_latitude_deg
%   (orbit never reaches the coverage latitude).

    if inclination_deg <= min_latitude_deg
        optimal_planes         = Inf;
        optimal_sats_per_plane = 0;
        total_sats             = Inf;
        return;
    end

    % WGS84 ellipsoid parameters (km)
    a = 6378.137;
    b = 6356.7523142;

    Rs = a + orbit_height_km;

    lat_rad = deg2rad(min_latitude_deg);
    sin2lat = sin(lat_rad)^2;
    cos2lat = cos(lat_rad)^2;

    Re_lat = sqrt((a^4 * cos2lat + b^4 * sin2lat) / ...
                  (a^2 * cos2lat + b^2 * sin2lat));

    % lambda_max: Earth Central Angle from nadir to edge of coverage footprint
    alpha    = asind((Re_lat / Rs) * cosd(min_elevation_deg));
    lmax_rad = deg2rad(180 - (90 + min_elevation_deg + alpha));

    % Longitude drift of counter-rotating track at critical latitude (Napier's Circle)
    delta_lon_rad = asin(tan(lat_rad) / tan(deg2rad(inclination_deg)));

    % Search over candidate values of satellites per plane
    min_s = ceil(pi / lmax_rad) + 1;
    max_s = 60;

    best_total             = Inf;
    optimal_sats_per_plane = 0;
    optimal_planes         = 0;

    for s = min_s:max_s
        S2 = pi / s;                    % half in-track spacing (ECA)

        if S2 >= lmax_rad
            continue;                   % street cannot be closed with this s
        end

        lambda_street = acos(cos(lmax_rad) / cos(S2));

        D_ctr  = 2 * lambda_street;         % max ECA budget across counter-rotating seam
        D_same = lambda_street + lmax_rad;  % max ECA budget between co-rotating planes

        % Convert ECA budgets to RAAN budgets on the sphere at latitude φ
        cos_arg_seam = (cos(D_ctr)  - sin2lat) / cos2lat;
        cos_arg_same = (cos(D_same) - sin2lat) / cos2lat;

        if cos_arg_seam < -1 || cos_arg_seam > 1 || ...
           cos_arg_same < -1 || cos_arg_same > 1
            continue;
        end

        R_seam_max = acos(cos_arg_seam) - 2 * delta_lon_rad;
        R_co_max   = acos(cos_arg_same);

        if R_seam_max <= 0
            continue;   % inclination drift exceeds available seam RAAN budget
        end

        P     = ceil((pi - R_seam_max) / R_co_max) + 1;
        total = P * s;

        if total < best_total
            best_total             = total;
            optimal_sats_per_plane = s;
            optimal_planes         = P;
        end
    end

    total_sats = best_total;
end
