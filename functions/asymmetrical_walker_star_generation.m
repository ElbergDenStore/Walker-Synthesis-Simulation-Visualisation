function sats = asymmetrical_walker_star_generation(orbit_height, inclination, planes, sats_per_plane)
        r_earth = 6378.14e3;
        a = r_earth + orbit_height;
        e = 0; % Circular orbit
        argPer = 0; 
        
        co_rotating_spacing = 180 / (planes - 1/3); % based on seamgap = 2/3 of normal gap

        
        % Calculate the Walker Phasing offset
        phase_shift = (Cfg.Phasing * 360) / Cfg.Total_sats;
        
        sat_array = []; % Temporary array to hold satellites
        
        for p = 1:planes
            % Calculate asymmetrical RAAN
            raan = (p - 1) * co_rotating_spacing;
            
            for s = 1:sats_per_plane
                % Calculate True Anomaly (position in the orbit + staggering)
                nu = (s - 1) * (360 / S) + (p - 1) * phase_shift;
                nu = mod(nu, 360); % Keep it cleanly within 0-360 degrees
                
                % Must match your existing "S4D_X" naming convention
                sat_num = (p - 1) * S + s;
                name = sprintf("S4D_%d", sat_num);
                
                % Create satellite and add to scenario
                new_sat = satellite(sc, a, e, inclination, raan, argPer, nu, ...
                    "Name", name, "OrbitPropagator", "sgp4");
                
                sat_array = [sat_array, new_sat]; 
            end
        end
        sats = sat_array; % Assign array to the expected variable