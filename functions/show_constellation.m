function show_constellation(Cfg, show_interactive, save_fig, out_dir, show_details)
    % Set default behaviors if you don't provide all inputs
    if nargin < 1
        Cfg = get_cfg(1000,"walkerdelta","big","short")
    end
    if nargin < 2, show_interactive = true; end
    if nargin < 3, save_fig = false; end
    if nargin < 4, out_dir = "figures"; end % Default to current folder
    if nargin < 5, show_details = false; end % Default to current folder


    sc = satelliteScenario;
    % EpochTime: if provided, satellites are initialised at this epoch so that their
    % orbital elements match the original simulation (not the viewer window start).
    % The viewer then jumps to Cfg.StartTime after opening.
    if isfield(Cfg, 'EpochTime')
        sc.StartTime = Cfg.EpochTime;
    else
        sc.StartTime = Cfg.StartTime;
    end
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;
    
    simTimes = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
    simTimes.TimeZone = 'UTC';
    r_earth = 6378.14e3;
    
    if Cfg.WalkerStar == true
        if isfield(Cfg, 'Lat_range_deg')
            min_lat_cov = min(Cfg.Lat_range_deg);
        else
            min_lat_cov = 0;
        end
        sats = asymmetrical_walker_star_generation(sc, Cfg.Orbit_height,  Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Min_elevation_UE, "two-body-keplerian", min_lat_cov);
    else
        sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, ...
        Cfg.Total_sats, ...
        Cfg.Num_planes, ...
        Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");
    end
    
    %% UEs Array
    if ~isfield(Cfg, 'Flat_UE_array') || ~isfield(Cfg.Flat_UE_array, 'Lats') || ~isfield(Cfg.Flat_UE_array, 'Lons')
        error('Cfg.Flat_UE_array with fields Lats and Lons is required. Generate UEs before calling coverage_simulator_function.');
    end

    UE_lats = Cfg.Flat_UE_array.Lats;
    UE_lons = Cfg.Flat_UE_array.Lons;
    %% Create the UEs Array
    NumUEs = length(UE_lats);
    % UEs = cell(NumUEs, 1);
    
    % for idx = 1:NumUEs
    %     UEs{idx}.Lat = UE_lats(idx);
    %     UEs{idx}.Lon = UE_lons(idx);
    %     UEs{idx}.Name = sprintf('UE%d', idx);
        
    %     % Add to scenario
    %     groundStation(sc, UEs{idx}.Lat, UEs{idx}.Lon, ...
    %         'Name', UEs{idx}.Name, 'MinElevationAngle', Cfg.Min_elevation_UE);
    % end
    %% Vectorized UE Array Creation
    % NumUEs = length(UE_lats);

    % 1. Generate all names at once as a string array (e.g., ["UE1", "UE2", ...])
    % ue_names = compose('UE%d', 1:NumUEs); 

    % 2. Create ALL ground stations in one single call
    % By passing arrays for lat/lon/names, MATLAB handles the loop internally in C++
    if ~isempty(UE_lats)
        groundStation(sc, UE_lats(:), UE_lons(:), ...
            'MinElevationAngle', Cfg.Min_elevation_UE);
    end


    if show_interactive || save_fig
        
        % 1. Launch the viewer FIRST so it's ready to receive graphics
        v = satelliteScenarioViewer(sc, 'ShowDetails', show_details);
        % Jump to the requested start time when epoch differs (e.g. diagnostic mode)
        if isfield(Cfg, 'EpochTime')
            v.CurrentTime = Cfg.StartTime;
        end

        for idx = 1:length(sc.GroundStations)
            sc.GroundStations(idx).ShowLabel = false; % even if showdetails is true, remove the UE labels
        end
        
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