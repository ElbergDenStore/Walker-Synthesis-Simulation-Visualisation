function sats = asymmetrical_walker_star_generation(sc, orbit_height, inclination, planes, sats_per_plane, min_elevation_deg)
    arguments
        sc (1,1) satelliteScenario
        orbit_height (1,1) double % Assumed to be in meters based on r_earth = 6378.14e3
        inclination (1,1) double
        planes (1,1) double {mustBeInteger, mustBePositive}
        sats_per_plane (1,1) double {mustBeInteger, mustBePositive}
        min_elevation_deg (1,1) double 
    end

    r_earth = 6378.14e3;
    a = r_earth + orbit_height;
    e = 0;
    argPer = 0;

    %  Calculate EXACT Seam Ratio for equal overlap ---
    orbit_height_km = orbit_height / 1000;
    Re_km = 6378.14;
    Rs_km = Re_km + orbit_height_km;

    alpha = asind((Re_km / Rs_km) * cosd(min_elevation_deg));
    lambda_max = deg2rad(180 - (90 + min_elevation_deg + alpha));
    
    S = (2 * pi) / sats_per_plane;
    lambda_street = acos(min(1, cos(lambda_max) / cos(S / 2)));
    
    D_maxCounter = 2 * lambda_street;
    D_maxSame = lambda_street + lambda_max;
    
    seam_ratio = D_maxCounter / D_maxSame; 

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
                'Name', name, 'OrbitPropagator', 'sgp4');

            sat_array = [sat_array, new_sat];
        end
    end

    sats = sat_array;
end