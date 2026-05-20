function [optimal_planes, optimal_sats_per_plane, total_sats] = calculate_walker_star(orbit_height_km, min_latitude_deg, min_elevation_deg)
    % WGS84 ellipsoid parameters (km)
    a = 6378.137;       % semi-major axis (equatorial radius)
    b = 6356.7523142;   % semi-minor axis (polar radius)

    % Orbital radius: circular orbit defined from the equatorial surface
    Rs = a + orbit_height_km;

    % WGS84 Earth radius at the target latitude
    lat_rad = deg2rad(min_latitude_deg);
    Re_lat = sqrt((a^4 * cos(lat_rad)^2 + b^4 * sin(lat_rad)^2) / ...
                  (a^2 * cos(lat_rad)^2 + b^2 * sin(lat_rad)^2));

    % Calculate Earth Central Angle (lambda_max) in radians
    alpha = asind((Re_lat / Rs) * cosd(min_elevation_deg));
    lambda_max = deg2rad(180 - (90 + min_elevation_deg + alpha));

    % Define search bounds for satellites per plane
    % Must have enough sats so S/2 is less than lambda_max to close the street
    min_sats_per_plane = ceil(pi / lambda_max) + 1; 
    max_sats_per_plane = 60; % Upper search limit

    best_total_sats = inf;
    optimal_sats_per_plane = 0;
    optimal_planes = 0;

    % Optimize by testing all viable integer values of sats per plane
    for s = min_sats_per_plane:max_sats_per_plane
        S = (2 * pi) / s;
        
        % Check if the street is physically viable (cosine ratio must be <= 1)
        if (S / 2) >= lambda_max
            continue; 
        end
        
        lambda_street = acos(cos(lambda_max) / cos(S / 2));
        
        D_maxCounter = 2 * lambda_street;
        D_maxSame = lambda_street + lambda_max;
        
        P = ceil((((a / Re_lat) * cosd(min_latitude_deg)*pi) - D_maxCounter) / D_maxSame) + 1;
        total_sats = P * s;
        
        % Store the configuration if it is the new minimum
        if total_sats < best_total_sats
            best_total_sats = total_sats;
            optimal_sats_per_plane = s;
            optimal_planes = P;
        end
    end

    total_sats = best_total_sats;
end