function sats = asymmetrical_walker_star_generation(sc, orbit_height, inclination, planes, sats_per_plane)
    arguments
        sc (1,1) satelliteScenario
        orbit_height (1,1) double
        inclination (1,1) double
        planes (1,1) double {mustBeInteger, mustBePositive}
        sats_per_plane (1,1) double {mustBeInteger, mustBePositive}
    end

    r_earth = 6378.14e3;
    a = r_earth + orbit_height;
    e = 0;
    argPer = 0;

    % Use the asymmetrical Walker-star seam spacing from the original design. % TODO - FIX THIS
    co_rotating_spacing = 180 / (planes - 1/3); % This is wrong and retarded
    in_plane_spacing = 360 / sats_per_plane; % correct
    total_sats = planes * sats_per_plane; % unused
    phase_shift = 180 / sats_per_plane; % is this correct?

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
