function sats = asymmetrical_walker_star_generation(sc, orbit_height, inclination, planes, sats_per_plane, min_elevation_deg, propagator, min_latitude_deg)
    arguments
        sc (1,1) satelliteScenario
        orbit_height (1,1) double % Assumed to be in meters based on r_earth = 6378.14e3
        inclination (1,1) double
        planes (1,1) double {mustBeInteger, mustBePositive}
        sats_per_plane (1,1) double {mustBeInteger, mustBePositive}
        min_elevation_deg (1,1) double
        propagator (1,1) string = "sgp4"
        min_latitude_deg (1,1) double = 0  % coverage reference latitude (deg); 0 = equatorial
    end

    r_earth = 6378.14e3;
    a = r_earth + orbit_height;

    % No eccentricity prewarping: SGP4 interprets the elements as Brouwer
    % mean elements, so e=0 already maps to a near-circular mean orbit.
    % The numerical propagator treats them as osculating, but prewarping does
    % not help there either.  The custom RK4+J2 propagator handles its own
    % frozen-orbit initialization internally.
    e      = 0;
    argPer = 0;

    %  Calculate EXACT Seam Ratio for equal overlap ---
    orbit_height_km = orbit_height / 1000;
    Re_eq_km = 6378.14;                             % equatorial radius for orbital mechanics
    Rs_km    = Re_eq_km + orbit_height_km;          % orbital radius from Earth centre
    % WGS84 surface radius at the minimum coverage latitude (for seam geometry)
    a_wgs = 6378.137; b_wgs = 6356.7523142;
    lat_r = deg2rad(min_latitude_deg);
    Re_km = sqrt((a_wgs^4*cos(lat_r)^2 + b_wgs^4*sin(lat_r)^2) / ...
                 (a_wgs^2*cos(lat_r)^2 + b_wgs^2*sin(lat_r)^2));

    alpha = asind((Re_km / Rs_km) * cosd(min_elevation_deg));
    lambda_max = deg2rad(180 - (90 + min_elevation_deg + alpha));
    
    S = (2 * pi) / sats_per_plane;
    lambda_street = acos(min(1, cos(lambda_max) / cos(S / 2)));
    
    D_maxCounter = 2 * lambda_street;
    D_maxSame = lambda_street + lambda_max;

    % Compute seam ratio in RAAN space (spherical law of cosines at coverage lat).
    % Using the ECA ratio D_ctr/D_same would slightly over-allocate RAAN to the
    % seam at high latitudes, producing a microscopic gap at boundary constellations.
    sin2lat = sin(lat_r)^2;
    cos2lat = cos(lat_r)^2;
    cos_seam = (cos(D_maxCounter) - sin2lat) / cos2lat;
    cos_co   = (cos(D_maxSame)   - sin2lat) / cos2lat;
    if cos_seam >= -1 && cos_seam <= 1 && cos_co > -1 && cos_co <= 1
        seam_ratio = acos(cos_seam) / acos(cos_co);   % exact RAAN-budget ratio
    else
        seam_ratio = D_maxCounter / D_maxSame;         % fallback: equatorial / low-lat
    end

    % --- Apply Spacing ---
    co_rotating_spacing = 180 / (planes - 1 + seam_ratio); 
    
    in_plane_spacing = 360 / sats_per_plane; 
    phase_shift = in_plane_spacing / 2; % Optimal staggered "brick wall" street

    sat_array = [];

    for p = 1:planes
        raan = (p - 1) * co_rotating_spacing;

        for s = 1:sats_per_plane
            nu = mod((s - 1) * in_plane_spacing + (p - 1) * phase_shift, 360);
            sat_num = (p - 1) * sats_per_plane + s;
            name = sprintf('S4D_%d', sat_num);

            new_sat = satellite(sc, a, e, inclination, raan, argPer, nu, ...
                'Name', name, 'OrbitPropagator', propagator);

            sat_array = [sat_array, new_sat];
        end
    end

    sats = sat_array;
end