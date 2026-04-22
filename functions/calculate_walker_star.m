function [optimal_planes, optimal_sats_per_plane, total_sats] = calculate_walker_star(orbit_height_km, min_latitude_deg, min_elevation_deg)
    Re = 6378.14;
    Rs = Re + orbit_height_km;

    % Calculate Earth Central Angle (lambda_max) in radians
    alpha = asind((Re / Rs) * cosd(min_elevation_deg));
    lambda_max = deg2rad(180 - (90 + min_elevation_deg + alpha));

    % Define search bounds for satellites per plane
    % Must have enough sats so S/2 is less than lambda_max to close the street
    min_sats_per_plane = ceil(pi / lambda_max) + 1; 
    max_sats_per_plane = 40; % Upper search limit

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
        
        P = ceil((((cosd(min_latitude_deg)*pi) - D_maxCounter) / D_maxSame) + 1);
        total_sats = P * s;
        
        % Store the configuration if it is the new minimum
        if total_sats < best_total_sats
            best_total_sats = total_sats;
            optimal_sats_per_plane = s;
            optimal_planes = P;
        end
    end
end