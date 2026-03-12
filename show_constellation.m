function show_constellation(Cfg, show_interactive, save_fig, out_dir)
    % 1. Set default behaviors if you don't provide all inputs
    if nargin < 1, 
        Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
        Cfg.StopTime   = datetime('1-Jun-2025 12:59:59', 'TimeZone', 'UTC');
        Cfg.SampleTime = 60; % seconds
        Cfg.Lat_vec = linspace(55, 85, 5);  
        Cfg.Lon_vec = linspace(-60, 30, 5);
        Cfg.Min_elevation_UE = 20;

        Cfg.Orbit_height = 1200e3;
        Cfg.Num_planes   = 4;
        Cfg.Inclination  = 75;
        Cfg.Total_sats   = 60;
        Cfg.Phasing      = 3; 
        Cfg.WalkerStar   = false; 
    end
    if nargin < 2, show_interactive = false; end
    if nargin < 3, save_fig = true; end
    if nargin < 4, out_dir = pwd; end % Default to current folder

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
        
        phase_shift = (Cfg.Phasing * 360) / Cfg.Total_sats;
        sat_array = []; 
        
        for p = 1:P
            raan = (p - 1) * corotating_gap;
            for s = 1:S
                nu = (s - 1) * (360 / S) + (p - 1) * phase_shift;
                nu = mod(nu, 360); 
                sat_num = (p - 1) * S + s;
                name = sprintf("S4D_%d", sat_num);
                
                new_sat = satellite(sc, a, e, inc, raan, argPer, nu, ...
                    "Name", name, "OrbitPropagator", "sgp4");
                sat_array = [sat_array, new_sat]; 
            end
        end
        sat = sat_array; 
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
   
    for idx = 1:NumUEs
        current_UE = UEs{idx};
        ue = groundStation(sc, current_UE.Lat, current_UE.Lon, ...
            'Name', current_UE.Name, 'MinElevationAngle', Cfg.Min_elevation_UE);
    end
    
    %% --- SMART VISUALIZATION LOGIC ---
    
    % Check if we are in a headless environment (Linux server w/ no display)
    % is_headless = ~usejava('desktop');
    % 
    % if is_headless
    %     fprintf('\n[i] Headless mode detected. 3D Globe rendering is skipped to prevent crashes.\n');
    %     fprintf('[i] Simulation math and scenario setup completed successfully.\n');
    %     return; % Exit the function safely here!
    % end
    
    % If we are here, we are on a computer with a GUI (like your Windows laptop)
    if show_interactive || save_fig
        
        % 1. Launch the viewer FIRST so it's ready to receive graphics
        v = satelliteScenarioViewer(sc, 'ShowDetails', false);
        
        % 2. Calculate and apply the sensors
        a = r_earth + Cfg.Orbit_height;
        A = Cfg.Min_elevation_UE + 90;
        b = r_earth;
        B = asin((b*sind(A))/a);
        max_view_angle = rad2deg(B)*2;
        
        sensors = conicalSensor(sat, 'MaxViewAngle', max_view_angle); 
        
        % CRITICAL FIX: Assign to a variable 'fov' so MATLAB doesn't delete it!
        fov = fieldOfView(sensors);
        % Explicitly draw the 3D orbital rings in space
        orb = orbit(sat);
        
        % 3. Position the camera
        target_lat = 57;
        target_lon = -9;
        target_alt = (r_earth + Cfg.Orbit_height) * 3;
        campos(v, target_lat, target_lon, target_alt);
        
        % 4. Force the GPU to draw the cones by nudging the time forward by 1 second
        % sc.SimulationTime = sc.StartTime + seconds(1); % throws error
        drawnow;
        if save_fig
            fprintf('Rendering 3D Globe... Please wait 8 seconds.\n');
            pause(8); % Give your GPU time to draw the sensor cones
            
            try
                % Grab the UI window by looking at ALL children, even hidden ones!
                viewer_fig = findall(allchild(0), 'Name', 'Satellite Scenario Viewer'); 
                
                % Fallback just in case the name varies by MATLAB version
                if isempty(viewer_fig)
                    viewer_fig = findall(allchild(0), 'Type', 'uifigure');
                end
                
                if ~isempty(viewer_fig)
                    % Make sure the output directory exists
                    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
                    
                    filename = fullfile(out_dir, sprintf('Constellation_3D_%d_Planes.png', Cfg.Num_planes));
                    
                    % Try exportapp (Standard for UIFigures)
                    exportapp(viewer_fig(1), filename);
                    fprintf('--> SUCCESS! Saved image to: %s\n', filename);
                    
                    % If you only wanted to save the image (like in a loop), close it automatically
                    if ~show_interactive
                        close(viewer_fig(1));
                    end
                else
                    fprintf('[!] Could not locate the hidden 3D viewer window to save.\n');
                end
            catch ME
                fprintf('[!] Failed to save image.\n');
                disp(ME.message);
            end
        end
    end
end