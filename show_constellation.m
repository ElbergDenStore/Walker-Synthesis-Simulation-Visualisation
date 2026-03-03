function show_constellation(Cfg)
    sc = satelliteScenario;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;
    
    simTimes = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
    simTimes.TimeZone = 'UTC';

    r_earth = 6378.14e3;
    
    if Cfg.WalkerStar == true
        sat = walkerStar(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, ...
        Cfg.Total_sats, ...
        Cfg.Num_planes, ...
        Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");
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