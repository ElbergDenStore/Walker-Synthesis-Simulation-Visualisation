function [planes, sats_per_plane, total_sats] = calculate_walker_star(orbit_height_km, min_latitude_deg, min_elevation_deg)
    Re = 6378.14; 
    Rs = Re + orbit_height_km; 
    
    % Circumference of the Earth at the minimum target latitude
    earth_O_at_lat = (cosd(min_latitude_deg) * Re) * 2 * pi;
    
    % Calculate Earth Central Angle (ECA) based on minimum allowable elevation
    alpha = asind((Re / Rs) * cosd(min_elevation_deg));
    ECA = deg2rad(180 - (90 + min_elevation_deg + alpha));
    
    % Hexagonal footprint dimensions
    hex_side = ECA * Re;
    max_normal_gap = 1.5 * hex_side;
    max_seam_gap = 2 * (hex_side/2);
    
    % planes = (earth circumference @ 55 lat / 2 ) / Beam coverage)
    % As the seam is smaller, ”planes – 1” needs to cover (earthcircumference @ 55 lat / 2) – seam coverage
    % Planes – 1 = ((earth circumference @ 55 lat / 2) – seam coverage) / beam coverage 
    planes = ceil((((earth_O_at_lat / 2) - max_seam_gap) / max_normal_gap) + 1);

    % Calculate Satellites per Plane
    sats_per_plane = ceil((2 * pi) / (sqrt(3) * ECA)); 
    
    % Output the total
    total_sats = planes * sats_per_plane;
end