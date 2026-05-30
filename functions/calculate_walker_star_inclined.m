function [optimal_planes, optimal_sats_per_plane, total_sats] = ...
    calculate_walker_star_inclined(orbit_height_km, min_latitude_deg, ...
                                   min_elevation_deg, inclination_deg)
% CALCULATE_WALKER_STAR_INCLINED  Minimum-satellite Walker Star for inclination < 90°.
%
%   Uses the spherical street-of-coverage model (Rider/Ballard/Walker) extended
%   for inclined orbits.  Uses a RAAN-budget formulation that matches the
%   Walker–SMAD formula in calculate_walker_star.m at i = 90°.
%
%   Physics:
%     At latitude φ a satellite track drifts in longitude (from Napier's Circle):
%         Δλ(φ,i) = arcsin( tan(φ) / tan(i) )
%     Counter-rotating planes at the seam share an ascending-node RAAN.  Their
%     tracks bow outward by ±Δλ, widening the coverage gap at the seam.
%
%     ECA budgets (D_ctr, D_same) are converted to RAAN longitude via the
%     WGS84-corrected latitude-circle approximation:
%         R = D · Re_lat / (a · cos φ)
%     This matches the SMAD total-width formula and avoids the optimism of
%     the exact great-circle inverse, which maps each ECA step to a slightly
%     larger longitude interval and produces analytically-marginal
%     configurations that fail in discrete-time simulation.
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

        % Convert ECA budgets to RAAN budgets using the WGS84-corrected
        % latitude-circle arc formula: R = D * Re_lat / (a * cos(φ)).
        %
        % The exact great-circle inverse (acos((cos(D)-sin²φ)/cos²φ)) maps
        % each ECA step to a slightly LARGER longitude interval than the
        % latitude-circle approximation, because the great-circle chord
        % between two equi-latitude points is shorter than the arc along
        % the parallel.  The surplus makes the formula optimistic: at large
        % footprints (low altitudes) it finds constellations whose coverage
        % margins are only a few tenths of a degree — analytically valid but
        % reliably missed by discrete-time simulation (660 s sample time).
        %
        % The latitude-circle formula matches the Walker–SMAD derivation in
        % calculate_walker_star.m (total width (a/Re_lat)·cos(φ)·π divided
        % by D_same ECA steps) and reproduces the numerically verified counts.
        lon_per_ECA = Re_lat / (a * cos(lat_rad));   % rad RAAN per rad ECA

        R_seam_max = D_ctr  * lon_per_ECA - 2 * delta_lon_rad;
        R_co_max   = D_same * lon_per_ECA;

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
