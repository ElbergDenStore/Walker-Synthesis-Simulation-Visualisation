function [optimal_planes, optimal_sats_per_plane, total_sats] = ...
    calculate_walker_star(orbit_height_km, min_latitude_deg, ...
                          min_elevation_deg, inclination_deg)
% CALCULATE_WALKER_STAR  Minimum-satellite Walker Star via street-of-coverage.
%
%   [optimal_planes, optimal_sats_per_plane, total_sats] = ...
%       calculate_walker_star(orbit_height_km, min_latitude_deg, ...
%                             min_elevation_deg, inclination_deg)
%
%   Closed-form (Rider/Ballard/Walker SMAD) sizing of the smallest Walker Star
%   that guarantees continuous coverage above min_latitude_deg at
%   min_elevation_deg.  Handles both polar (i = 90 deg) and inclined orbits.
%
%   Inputs
%     orbit_height_km    Orbital altitude above the equatorial surface (km)
%     min_latitude_deg   Lowest latitude that must be served (deg)
%     min_elevation_deg  Minimum service elevation angle (deg)
%     inclination_deg    (optional) Orbital inclination (deg). Default 90.
%
%   Outputs
%     optimal_planes         Number of orbital planes (P)
%     optimal_sats_per_plane Satellites per plane (S)
%     total_sats             Minimum total satellite count (P*S)
%
%   Physics:
%     At latitude phi a satellite track drifts in longitude (Napier's Circle):
%         delta_lon(phi, i) = arcsin( tan(phi) / tan(i) )
%     Counter-rotating planes at the seam share an ascending-node RAAN; their
%     tracks bow outward by +/- delta_lon, widening the coverage gap at the
%     seam.  ECA budgets (D_ctr, D_same) are converted to RAAN longitude via
%     the WGS84-corrected latitude-circle approximation:
%         R = D * Re_lat / (a * cos(phi))
%     which matches the SMAD total-width formula and avoids the optimism of
%     the exact great-circle inverse (which yields analytically-marginal
%     configurations that fail in discrete-time simulation).
%
%   At i = 90 deg the drift term delta_lon vanishes and the RAAN-budget
%   formulation reduces exactly to the classic polar SMAD street-of-coverage
%   formula, so 3-argument calls reproduce the polar result bit-for-bit.
%
%   Returns [Inf, 0, Inf] if inclination_deg <= min_latitude_deg
%   (the orbit never reaches the required coverage latitude).

    if nargin < 4 || isempty(inclination_deg)
        inclination_deg = 90;   % default: polar Walker Star
    end

    if inclination_deg <= min_latitude_deg
        optimal_planes         = Inf;
        optimal_sats_per_plane = 0;
        total_sats             = Inf;
        return;
    end

    % WGS84 ellipsoid parameters (km)
    a = 6378.137;       % semi-major axis (equatorial radius)
    b = 6356.7523142;   % semi-minor axis (polar radius)

    % Orbital radius: circular orbit defined from the equatorial surface
    Rs = a + orbit_height_km;

    % WGS84 Earth radius at the target latitude
    lat_rad = deg2rad(min_latitude_deg);
    sin2lat = sin(lat_rad)^2;
    cos2lat = cos(lat_rad)^2;
    Re_lat  = sqrt((a^4 * cos2lat + b^4 * sin2lat) / ...
                   (a^2 * cos2lat + b^2 * sin2lat));

    % Earth Central Angle to the edge of the coverage footprint (rad)
    alpha      = asind((Re_lat / Rs) * cosd(min_elevation_deg));
    lambda_max = deg2rad(180 - (90 + min_elevation_deg + alpha));

    % Longitude drift of a counter-rotating track at the critical latitude
    % (Napier's Circle). Zero for polar orbits (i = 90 deg).
    delta_lon_rad = asin(tan(lat_rad) / tan(deg2rad(inclination_deg)));

    % RAAN longitude per unit ECA (latitude-circle arc, WGS84-corrected).
    % At i = 90 this makes the formula below equal the polar SMAD form
    % P = ceil(((a/Re_lat)*cos(phi)*pi - D_ctr)/D_same) + 1.
    lon_per_ECA = Re_lat / (a * cos(lat_rad));

    % Search bounds for satellites per plane: need enough sats so S/2 < lambda_max
    min_sats_per_plane = ceil(pi / lambda_max) + 1;
    max_sats_per_plane = 60;

    best_total_sats        = inf;
    optimal_sats_per_plane = 0;
    optimal_planes         = 0;

    % Optimize by testing all viable integer values of sats per plane
    for s = min_sats_per_plane:max_sats_per_plane
        S = (2 * pi) / s;

        % Street must be physically closable (cosine ratio must be <= 1)
        if (S / 2) >= lambda_max
            continue;
        end

        lambda_street = acos(cos(lambda_max) / cos(S / 2));

        D_maxCounter = 2 * lambda_street;          % ECA budget across counter-rotating seam
        D_maxSame    = lambda_street + lambda_max;  % ECA budget between co-rotating planes

        % Convert ECA budgets to RAAN budgets; subtract the inclination drift
        % from the seam budget (it widens the counter-rotating gap).
        R_seam_max = D_maxCounter * lon_per_ECA - 2 * delta_lon_rad;
        R_co_max   = D_maxSame    * lon_per_ECA;

        if R_seam_max <= 0
            continue;   % inclination drift exceeds available seam RAAN budget
        end

        P          = ceil((pi - R_seam_max) / R_co_max) + 1;
        total_sats = P * s;

        % Store the configuration if it is the new minimum
        if total_sats < best_total_sats
            best_total_sats        = total_sats;
            optimal_sats_per_plane = s;
            optimal_planes         = P;
        end
    end

    total_sats = best_total_sats;
end