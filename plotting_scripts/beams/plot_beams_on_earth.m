close all force;
clear variables;
clc;

%% Output Directory
out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'plotting_scripts/figures', 'beams_on_earth');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% 1. User Config
lat = 57;
lon = 9.3;
alt_m = 550e3;

f_hz = 12e9;
G_tx_dBi = 37;
min_elev_deg = 20;
frf = 3;

% --- Figure Style ---
fig_size   = [500, 500];   % [width, height] in pixels
font_size  = 14;           % axis labels, tick labels, title

%% 2. Scenario & Satellite
startTime = datetime('now');
stopTime = startTime + minutes(1);
sc = satelliteScenario(startTime, stopTime, 10);

posData = [lat, lon, alt_m; lat, lon, alt_m];
tt = timetable([startTime; stopTime], posData);
sat = satellite(sc, tt, 'Name', 'Multi-Beam LEO', 'CoordinateFrame', 'geographic');

%% 3. Use calculate_hexagonal_beams (single geometry source)
BeamGrid = calculate_hexagonal_beams(G_tx_dBi, f_hz, alt_m, min_elev_deg, [], frf);

b_u = BeamGrid.u_center(:);
b_v = BeamGrid.v_center(:);
num_beams = BeamGrid.num_beams;
uv_3dB = 2 * (1.391 / (BeamGrid.Nu * (pi/2)));
beamwidth_deg = 2 * asind(uv_3dB / 2);

% Convert uv to steering-space angles used by MATLAB sensor mounting.
[roll_deg, pitch_deg] = uv_to_matlab_mounting_angles(b_u, b_v);

% Viewer is circular by sensor model; we still size each beam with scan-loss stretch.
eta_deg = asind(sqrt(b_u.^2 + b_v.^2));
plot_width_deg = beamwidth_deg ./ cosd(eta_deg);

fprintf('Generating %d beams from calculate_hexagonal_beams...\n', num_beams);
names = "Beam_" + string(1:num_beams);
mounting_angles = [zeros(num_beams, 1), pitch_deg, roll_deg]';

% tic
% conicalSensor(sat, 'Name', names, ...
%                    'MaxViewAngle', plot_width_deg, ...
%                    'MountingAngles', mounting_angles);
% fprintf('Viewer sensors created in %.3f s\n', toc);

Re_km = 6371;
h_km = alt_m / 1000;
eta_max_deg = asind((Re_km / (Re_km + h_km)) * cosd(min_elev_deg));
% conicalSensor(sat, 'Name', "Min_Elevation_UE", 'MaxViewAngle', eta_max_deg * 2);

% %% 4. Viewer output - super cool 3d interactive plot, but number 6 is more nice
% v = satelliteScenarioViewer(sc, 'ShowDetails', false);
% fieldOfView(sat.ConicalSensors);
% campos(v, lat, lon, alt_m * 3);

% % Zoom pass: show only the beams close to nadir and save a clean snapshot.
% nadir_cutoff_deg = max(beamwidth_deg, 8);
% nadir_mask = eta_deg <= nadir_cutoff_deg;
% if any(nadir_mask)
%     fieldOfView(sat.ConicalSensors(nadir_mask));
%     campos(v, lat, lon, max(alt_m * 0.8, 50e3));
%     drawnow;
%     pause(2);
%     save_viewer_snapshot('screenshots/beams_viewer_nadir_zoom.png');
%     fprintf('Saved viewer snapshot: screenshots/beams_viewer_nadir_zoom.png\n');
% else
%     fprintf('No nadir beams matched the cutoff %.1f deg.\n', nadir_cutoff_deg);
% end

%% 5. Ellipsoidal plot in steering-angle space (exported)
theta = linspace(0, 2*pi, 80);
b_eta = eta_deg;
b_az = atan2d(b_v, b_u);
b_x = b_eta .* cosd(b_az);
b_y = b_eta .* sind(b_az);

fig = figure('Color', 'w', 'Name', 'Ellipsoidal Beam Footprint (Steering Space)', ...
             'Position', [100, 100, fig_size(1), fig_size(2)]);
hold on; axis equal; box on; grid on;

for b = 1:num_beams
    eta_b = b_eta(b);
    az_b = b_az(b);

    if eta_b < 1e-3
        map_dist = 1;
    else
        map_dist = (eta_b * pi/180) / sind(eta_b);
    end

    r_major = (beamwidth_deg / 2) / cosd(eta_b);
    r_minor = (beamwidth_deg / 2) * map_dist;

    x_local = r_major * cos(theta);
    y_local = r_minor * sin(theta);

    x_rot = x_local * cosd(az_b) - y_local * sind(az_b);
    y_rot = x_local * sind(az_b) + y_local * cosd(az_b);

    plot(b_x(b) + x_rot, b_y(b) + y_rot, 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
end

plot(eta_max_deg * cos(theta), eta_max_deg * sin(theta), 'k--', 'LineWidth', 1.8);
xlabel('X Steering Angle (deg off nadir)', 'FontSize', font_size);
ylabel('Y Steering Angle (deg off nadir)', 'FontSize', font_size);
title({'Beam Footprints from calculate\_hexagonal\_beams', ...
       sprintf('f = %.1f GHz  |  Gtx = %.1f dBi  |  FRF = %d', f_hz/1e9, G_tx_dBi, frf)}, ...
    'FontSize', font_size);
set(gca, 'FontSize', font_size);

exportgraphics(fig, fullfile(out_dir, 'beams_ellipsoidal_steering_space.png'), 'Resolution', 600);
fprintf('Saved plot: %s\n', fullfile(out_dir, 'beams_ellipsoidal_steering_space.png'));

%% 6. Phased Array Projection on Geographic Map
% This projects the steering-space ellipses onto the WGS84 Earth surface

% Earth geometry
Re_km = 6371;
h_km = alt_m / 1000;
rsat_km = Re_km + h_km;

fig_map = figure('Color', 'w', 'Name', 'Geographic Phased Array Footprint', ...
                 'Position', [100, 100, fig_size(1), fig_size(2)]);
gx = geoaxes('Basemap', 'satellite', 'FontSize', font_size);
hold on;

% 6a. Plot the individual beams
theta = linspace(0, 2*pi, 80);
fprintf('Projecting %d beams onto the Earth surface...\n', num_beams);

for b = 1:num_beams
    eta_b = b_eta(b);
    az_b = b_az(b);
    
    % Scan loss calculation
    if eta_b < 1e-3
        map_dist = 1;
    else
        map_dist = (eta_b * pi/180) / sind(eta_b);
    end
    
    r_major = (beamwidth_deg / 2) / cosd(eta_b);
    r_minor = (beamwidth_deg / 2) * map_dist;
    
    x_local = r_major * cos(theta);
    y_local = r_minor * sin(theta);
    
    x_rot = x_local * cosd(az_b) - y_local * sind(az_b);
    y_rot = x_local * sind(az_b) + y_local * cosd(az_b);
    
    x_contour = b_x(b) + x_rot;
    y_contour = b_y(b) + y_rot;
    
    eta_contour = sqrt(x_contour.^2 + y_contour.^2);
    az_compass = atan2d(x_contour, y_contour); 
    
    sin_arg = (rsat_km / Re_km) .* sind(eta_contour);
    valid_pts = sin_arg <= 1; 
    
    lambda_deg = nan(size(eta_contour));
    lambda_deg(valid_pts) = asind(sin_arg(valid_pts)) - eta_contour(valid_pts);
    
    [lat_contour, lon_contour] = reckon(lat, lon, lambda_deg, az_compass);
    geoplot(gx, lat_contour, lon_contour, 'Color', [1 0.2 0.8], 'LineWidth', 1.0);
end

% 6b. Plot the 20-Degree Minimum Elevation Boundary
az_boundary = linspace(0, 360, 360);
sin_arg_bound = (rsat_km / Re_km) .* sind(eta_max_deg);
% Ensure the boundary doesn't exceed the Earth's horizon
if sin_arg_bound <= 1
    lambda_max_deg = asind(sin_arg_bound) - eta_max_deg;
    [lat_bound, lon_bound] = reckon(lat, lon, lambda_max_deg, az_boundary);
    geoplot(gx, lat_bound, lon_bound, 'Color', [1 0.9 0.2], 'LineStyle', '--', 'LineWidth', 2.5);
else
    % Fallback if 20-deg elevation somehow misses Earth
    lambda_max_deg = asind(1) - asind(Re_km/rsat_km); 
end

title(gx, {'Phased Array Footprints on Earth', ...
           sprintf('h = %d km  |  f = %.1f GHz  |  G = %.1f dBi', alt_m*1e-3, f_hz/1e9, G_tx_dBi)}, ...
    'FontSize', font_size);

% 6c. Export 3 specific zoom levels
fprintf('Saving figures (pausing between views to allow map tiles to load)...\n');

% Shot 1: Full Coverage Region
geolimits(gx, [lat - lambda_max_deg - 2, lat + lambda_max_deg + 2], ...
              [lon - lambda_max_deg*1.5 - 2, lon + lambda_max_deg*1.5 + 2]);
drawnow;
pause(5); % Wait 5 seconds for the massive full-res map to download
exportgraphics(fig_map, fullfile(out_dir, 'beams_01_full_coverage.png'), 'Resolution', 300);
fprintf('Saved: beams_01_full_coverage.png\n');

% Shot 2: Northern Jutland / Aalborg Zoom (~50 beams focus)
geolimits(gx, [lat - 0.8, lat + 0.8], [lon - 1.5, lon + 1.5]);
drawnow;
pause(3); % Wait 3 seconds for local zoom tiles
exportgraphics(fig_map, fullfile(out_dir, 'beams_02_jutland_zoom.png'), 'Resolution', 300);
fprintf('Saved: beams_02_jutland_zoom.png\n');

% Shot 3: Outskirts / Helsinki Zoom
helsinki_lat = 60.1695;
helsinki_lon = 24.9354;
% Framed to show Helsinki and the array edge dropping into the Gulf of Finland
geolimits(gx, [helsinki_lat - 1.2, helsinki_lat + 1.2], [helsinki_lon - 2.5, helsinki_lon + 2.5]);
drawnow;
pause(3); % Wait 3 seconds for Helsinki tiles
exportgraphics(fig_map, fullfile(out_dir, 'beams_03_helsinki_zoom.png'), 'Resolution', 300);
fprintf('Saved: beams_03_helsinki_zoom.png\n');

function [roll_deg, pitch_deg] = uv_to_matlab_mounting_angles(u, v)
% MATLAB applies mounting rotations in a sequence where roll influences
% pitch mapping in uv-space. Keep this exact order for correct pointing.

    v_clamped = min(max(v, -1), 1);
    roll_deg = asind(v_clamped);

    denom = cosd(roll_deg);
    denom = max(denom, 1e-8);
    pitch_arg = u ./ denom;
    pitch_arg = min(max(pitch_arg, -1), 1);
    pitch_deg = asind(pitch_arg);
end

function save_viewer_snapshot(out_file)
    viewer_fig = findall(allchild(0), 'Name', 'Satellite Scenario Viewer');
    if isempty(viewer_fig)
        viewer_fig = findall(allchild(0), 'Type', 'uifigure');
    end

    if isempty(viewer_fig)
        warning('Could not locate the satellite scenario viewer window for export.');
        return;
    end

    % Make the exported image sharper by enlarging the viewer before capture.
    try
        viewer_fig(1).WindowState = 'maximized';
    catch
        try
            viewer_fig(1).Position = [50 50 2200 1400];
        catch
        end
    end

    drawnow;
    pause(0.5);

    out_dir = fileparts(out_file);
    if ~isempty(out_dir) && ~isfolder(out_dir)
        mkdir(out_dir);
    end

    exportapp(viewer_fig(1), out_file);
end