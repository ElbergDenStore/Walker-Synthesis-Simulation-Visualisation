function sats = generate_walker_star_scenario(sc, orbit_height, inclination, planes, sats_per_plane, min_elevation_deg, propagator, min_latitude_deg)
% GENERATE_WALKER_STAR_SCENARIO  Add a Walker Star constellation to a scenario.
%   Populates the satelliteScenario SC with a Walker Star (asymmetric seam-ratio
%   RAAN spacing and brick-wall in-plane stagger) and returns the satellite
%   handles.  ORBIT_HEIGHT is in metres; INCLINATION, MIN_ELEVATION_DEG and
%   MIN_LATITUDE_DEG in degrees.  PROPAGATOR selects "sgp4" (default), etc.
    arguments
        sc (1,1) satelliteScenario
        orbit_height (1,1) double % Assumed to be in meters based on r_earth = 6378.137e3
        inclination (1,1) double
        planes (1,1) double {mustBeInteger, mustBePositive}
        sats_per_plane (1,1) double {mustBeInteger, mustBePositive}
        min_elevation_deg (1,1) double
        propagator (1,1) string = "sgp4"
        min_latitude_deg (1,1) double = 0  % coverage reference latitude (deg); 0 = equatorial
    end

    r_earth = 6378.137e3;
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
    Re_eq_km = 6378.137;                            % equatorial radius for orbital mechanics, WGS84
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
                'Name', name, 'OrbitPropagator', propagator);

            sat_array = [sat_array, new_sat];
        end
    end

    sats = sat_array;
end