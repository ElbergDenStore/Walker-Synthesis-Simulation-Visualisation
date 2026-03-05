function show_constellation(Cfg)
    sc = satelliteScenario;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;
    
    simTimes = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
    simTimes.TimeZone = 'UTC';

    r_earth = 6378.14e3;
    
    if Cfg.WalkerStar == true
        % CUSTOM ASYMMETRICAL WALKER STAR
        a = r_earth + Cfg.Orbit_height;
        e = 0; % Circular orbit
        inc = Cfg.Inclination;
        argPer = 0; 
        
        S = Cfg.Total_sats / Cfg.Num_planes; % Sats per plane
        P = Cfg.Num_planes; % Keep as total planes over 180 degrees
        

        corotating_gap = (180 - Cfg.Seam_gap) / (P - 1);
        
        % Calculate the Walker Phasing offset
        phase_shift = (Cfg.Phasing * 360) / Cfg.Total_sats;
        
        sat_array = []; % Temporary array to hold satellites
        
        for p = 1:P
            % Calculate asymmetrical RAAN
            raan = (p - 1) * corotating_gap;
            
            for s = 1:S
                % Calculate True Anomaly (position in the orbit + staggering)
                nu = (s - 1) * (360 / S) + (p - 1) * phase_shift;
                nu = mod(nu, 360); % Keep it cleanly within 0-360 degrees
                
                % Must match your existing "S4D_X" naming convention
                sat_num = (p - 1) * S + s;
                name = sprintf("S4D_%d", sat_num);
                
                % Create satellite and add to scenario
                new_sat = satellite(sc, a, e, inc, raan, argPer, nu, ...
                    "Name", name, "OrbitPropagator", "sgp4");
                
                sat_array = [sat_array, new_sat]; 
            end
        end
        sat = sat_array; % Assign array to the expected variable
    else
        sat = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, ...
        Cfg.Total_sats, ...
        Cfg.Num_planes, ...
        Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");
    end

    %% UE Grid Setup
    [LonGrid, LatGrid] = meshgrid(Cfg.Lon_vec, Cfg.Lat_vec);
    NumUEs = numel(LonGrid);
    UEs = cell(NumUEs, 1);
    for idx = 1:NumUEs
        UEs{idx}.Lat = LatGrid(idx);
        UEs{idx}.Lon = LonGrid(idx);
        UEs{idx}.Name = sprintf('UE%dLat%.0f', idx, UEs{idx}.Lat);
    end
   
    tic
    for idx = 1:NumUEs
        current_UE = UEs{idx};
        ue = groundStation(sc, current_UE.Lat, current_UE.Lon, ...
            'Name', current_UE.Name, 'MinElevationAngle', Cfg.Min_elevation_UE);
    end
    % a/Sin(A) = b/Sin(B)
    % sin(B) = b*sin(A)/a
    a = r_earth + Cfg.Orbit_height;
    A = Cfg.Min_elevation_UE + 90;
    b = r_earth;
    B = asin((b*sind(A))/a)
    max_view_angle = rad2deg(B)*2

    sensors = conicalSensor(sat, 'MaxViewAngle', max_view_angle); 
    fieldOfView(sensors);
    
    satelliteScenarioViewer(sc);
end