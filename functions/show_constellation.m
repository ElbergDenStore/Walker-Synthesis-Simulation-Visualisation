function show_constellation(Cfg, show_interactive, save_fig, out_dir)
    % Set default behaviors if you don't provide all inputs
    if nargin < 1
        Cfg.StartTime  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
        Cfg.StopTime   = datetime('1-Jun-2025 12:59:59', 'TimeZone', 'UTC');
        Cfg.SampleTime = 60; % seconds
        Cfg.Lat_vec = linspace(55, 85, 5);  
        Cfg.Lon_vec = linspace(-60, 30, 2);
        Cfg.Min_elevation_UE = 20;

        Cfg.Orbit_height = 1000e3;
        Cfg.Num_planes   = 5;
        Cfg.Sats_per_plane   = 13;
        Cfg.Inclination  = 87;
        Cfg.Total_sats   = Cfg.Num_planes * Cfg.Sats_per_plane;
        Cfg.WalkerStar     = true;
        Cfg.Phasing        = Cfg.Num_planes/2;
    end
    if nargin < 2, show_interactive = true; end
    if nargin < 3, save_fig = false; end
    if nargin < 4, out_dir = pwd; end % Default to current folder

    sc = satelliteScenario;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;
    
    simTimes = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
    simTimes.TimeZone = 'UTC';
    r_earth = 6378.14e3;
    
    if Cfg.WalkerStar == true
        sats = asymmetrical_walker_star_generation(orbit_height, inclination, planes, sats_per_plane);
    else
        sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, ...
        Cfg.Total_sats, ...
        Cfg.Num_planes, ...
        Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");
    end

    [UE_lats_flat, UE_lons_flat] = generate_equal_ish_area_UEs(Cfg.Lat_vec, Cfg.Lon_vec);
    
    %% Create the UEs Array
    NumUEs = length(UE_lats_flat);
    UEs = cell(NumUEs, 1);
    
    for idx = 1:NumUEs
        UEs{idx}.Lat = UE_lats_flat(idx);
        UEs{idx}.Lon = UE_lons_flat(idx);
        UEs{idx}.Name = sprintf('UE%d', idx);
        
        % Add to scenario
        groundStation(sc, UEs{idx}.Lat, UEs{idx}.Lon, ...
            'Name', UEs{idx}.Name, 'MinElevationAngle', Cfg.Min_elevation_UE);
    end

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
        
        sensors = conicalSensor(sats, 'MaxViewAngle', max_view_angle); 
        
        % Assign to a variable so matlab doesn't delete it
        fov = fieldOfView(sensors);
        % Explicitly draw the 3D orbital rings in space
        orb = orbit(sats);
        
        % 3. Position the camera
        target_lat = 57;
        target_lon = -9;
        % target_alt = (r_earth + Cfg.Orbit_height) * 2;
        target_alt = (r_earth + 1000) * 2; % Same camera height for all runs. Does not change much, but anyway
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