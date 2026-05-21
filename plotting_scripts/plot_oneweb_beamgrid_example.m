% --- 1. Setup the Environment ---
close all force;
clear variables;
clc;

path_setup()
% Call the function we created to get the beam definitions
BeamGrid = calculate_OneWeb_beams();

% --- 2. Create the Dense u, v Grid ---
% u is cross-track (wide beam), v is along-track (narrow steering)
% v needs to cover at least +/- 0.45 to see the +/- 25 degree beams
u_vec = linspace(-0.9, 0.9, 1000); 
v_vec = linspace(-0.9, 0.9, 1000); 
[U, V] = meshgrid(u_vec, v_vec);

num_beams = length(BeamGrid.v_center);

% --- 3. Evaluate Array Factor for All Beams ---
% We will store the normalized gain of every beam at every pixel in a 3D matrix
Normalized_Gain_3D = zeros(size(U, 1), size(U, 2), num_beams);

% Element Factor (Angle to Nadir)
% theta is the angle from nadir: sin^2(theta) = u^2 + v^2
sin_theta_sq = U.^2 + V.^2;
cos_theta = zeros(size(U));
cos_theta(sin_theta_sq <= 1) = sqrt(1 - sin_theta_sq(sin_theta_sq <= 1));
% Compute element power roll-off based on Cos_exponent
EF_dB = 10 * log10(max(cos_theta.^BeamGrid.Cos_exponent, 1e-10));

for b = 1:num_beams
    % 1. Distance from beam center (Mechanical face)
    du = U - BeamGrid.u_center(b);
    dv = V - BeamGrid.v_center(b);
    
    % Prevent exact zero division
    du(du == 0) = eps; 
    dv(dv == 0) = eps;
    
    % 2. Calculate Standard Array Factor (Assumes 0.5 lambda spacing)
    Nu = BeamGrid.Nu;
    Nv = BeamGrid.Nv;
    
    AF_u = sin(Nu * (pi/2) .* du) ./ (Nu .* sin((pi/2) .* du));
    AF_v = sin(Nv * (pi/2) .* dv) ./ (Nv .* sin((pi/2) .* dv));
    
    % Array Factor in dB
    AF_dB = 20 * log10(abs(AF_u .* AF_v));
    
    % 3. ELEMENT FACTOR (Centered on the beam, NOT Nadir)
    % Because the stick is mechanically tilted, theta is just the distance from beam center
    rho_sq = du.^2 + dv.^2;
    rho_sq(rho_sq > 1) = 1; 
    
    cos_theta = sqrt(1 - rho_sq);
    EF_loss_dB = 10 * BeamGrid.Cos_exponent * log10(max(cos_theta, eps));
    
    % 4. Total Normalized Gain
    % Because the element is broadside to the beam, the peak is naturally exactly 0 dB.
    loss_dB = AF_dB + EF_loss_dB;
    loss_dB = max(loss_dB, -60); % Cap at -60dB for visual clarity
    
    % Store normalized gain
    Normalized_Gain_3D(:, :, b) = loss_dB;
end

% --- 4. Calculate Serving Beam and Interference ---
% The serving beam at any pixel is the one providing the maximum normalized gain
[Max_Normalized_Map, Serving_Idx_Map] = max(Normalized_Gain_3D, [], 3);

% Calculate neighbor interference
Interference_Linear = zeros(size(U));

for b = 1:num_beams
    % Find all pixels where beam 'b' is the server
    server_mask = (Serving_Idx_Map == b);
    
    if ~any(server_mask, 'all')
        continue;
    end
    
    % Get valid neighbors for this beam
    neighbors = BeamGrid.neighbor_idx(b, :);
    neighbors = neighbors(~isnan(neighbors));
    
    % Sum the linear power from all neighboring beams at those pixels
    for n = neighbors
        neighbor_gain_dB = Normalized_Gain_3D(:, :, n);
        Interference_Linear(server_mask) = Interference_Linear(server_mask) + 10.^(neighbor_gain_dB(server_mask) ./ 10);
    end
end

% Convert interference back to dB
Interference_Map_dB = 10 * log10(Interference_Linear);
Interference_Map_dB(Interference_Linear == 0) = -Inf; % Handle areas with 0 neighbors

% --- 5. Plotting ---
% v_margin = 0.05;
u_lims = [-0.6, 0.6]; % Fan beams are very wide in u (cross-track)
v_lims = [-0.6, 0.6];

% Ensure the figures directory exists
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures', 'beamgrid');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Plot 1: 3D Surf Plot
fig1 = figure('Name', 'OneWeb Normalized Gain 3D Surf', 'Position', [100, 100, 600, 400]);
surf(U, V, Max_Normalized_Map, 'EdgeColor', 'none');
colormap(gca, 'turbo');
caxis([-15, 0]); % Bound the colormap limits for readability
colorbar;
hold on;

% Plot the beam centers on the surf surface
scatter3(BeamGrid.u_center, BeamGrid.v_center, zeros(1, num_beams), 50, 'k', 'filled');

xlabel('u (Cross-Track Direction Cosine)');
ylabel('v (Along-Track Direction Cosine)');
zlabel('Normalized Gain (dB)');
title('3D Map of Maximum Normalized Gain');
axis tight;
xlim(u_lims);
ylim(v_lims);
zlim([-15, 0]);
view(-30, 45); % adjust angle for 3D

fig1_filename = fullfile(fig_dir, 'oneweb_surf_gain');
exportgraphics(fig1, [fig1_filename, '.png'], 'Resolution', 300);

%% Plot 2: 2D Contours Only
fig2 = figure('Name', 'OneWeb Beams -3dB Contours', 'Position', [100, 100, 600, 400]);
hold on;
for b = 1:num_beams
    contour(U, V, Normalized_Gain_3D(:,:,b), [-3, -3], 'LineWidth', 1.5, 'LineColor', 'b');
end
plot(BeamGrid.u_center, BeamGrid.v_center, 'k.', 'MarkerSize', 12);

xlabel('u (Cross-Track Direction Cosine)');
ylabel('v (Along-Track Direction Cosine)');
title('OneWeb Beams -3dB Contours');
axis equal; axis tight;
xlim(u_lims);
ylim(v_lims);
grid on;

fig2_filename = fullfile(fig_dir, 'oneweb_contours_only');
exportgraphics(fig2, [fig2_filename, '.png'], 'Resolution', 300);

%% Plot 3: Geographic Projection over Denmark
fig3 = figure('Color', 'w', 'Name', 'Geographic Projection Denmark', 'Position', [100, 100, 600, 400]);
gx = geoaxes('Basemap', 'satellite');
hold(gx, 'on');

% Satellite location parameters (match beams_on_earth limits style)
lat = 57;
lon = 9.3;
alt_m = 1200e3; % Approx OneWeb height
Re_km = 6371;
rsat_km = Re_km + (alt_m / 1000);

% Extract -3dB contours manually to project them
for b = 1:num_beams
    % Extract the contour from the matrix
    c_matrix = contourc(u_vec, v_vec, Normalized_Gain_3D(:,:,b), [-3, -3]);
    
    idx = 1;
    while idx < size(c_matrix, 2)
        n_pts = c_matrix(2, idx);
        
        c_u = c_matrix(1, idx+1 : idx+n_pts);
        c_v = c_matrix(2, idx+1 : idx+n_pts);
        
        % Convert u, v to steering angles
        eta_contour = asind(sqrt(c_u.^2 + c_v.^2));
        az_compass = atan2d(c_u, c_v); 
        
        sin_arg = (rsat_km / Re_km) .* sind(eta_contour);
        valid_pts = sin_arg <= 1;
        
        if any(valid_pts)
            lambda_deg = asind(sin_arg(valid_pts)) - eta_contour(valid_pts);
            [lat_contour, lon_contour] = reckon(lat, lon, lambda_deg, az_compass(valid_pts));
            geoplot(gx, lat_contour, lon_contour, 'Color', [1 0.2 0.8], 'LineWidth', 1.5);
        end
        
        idx = idx + n_pts + 1;
    end
end

% Project beam centers
eta_center = asind(sqrt(BeamGrid.u_center.^2 + BeamGrid.v_center.^2));
az_center = atan2d(BeamGrid.u_center, BeamGrid.v_center);
sin_arg_center = (rsat_km / Re_km) .* sind(eta_center);
valid_centers = sin_arg_center <= 1;
if any(valid_centers)
    lambda_center_deg = asind(sin_arg_center(valid_centers)) - eta_center(valid_centers);
    [lat_center, lon_center] = reckon(lat, lon, lambda_center_deg, az_center(valid_centers));
    geoplot(gx, lat_center, lon_center, 'k.', 'MarkerSize', 12);
end

% Wait for tiles to load
fprintf('Saving geographic screenshots (pausing for map tiles to load)...\n');

% Shot 1: Full Europe / Wide Coverage
geolimits(gx, [46, 68], [-10, 28]);
% geolimits(gx, [-20, 20], [-20, 20]);
title(gx, 'OneWeb Beams -3dB footprints');
drawnow;
pause(5); % wait for wide view tiles
fig3_filename_wide = fullfile(fig_dir, 'oneweb_projection_wide');
exportgraphics(fig3, [fig3_filename_wide, '.png'], 'Resolution', 300);

% Shot 2: Northern Jutland / Aalborg Zoom
geolimits(gx, [lat - 0.8, lat + 0.8], [lon - 1.5, lon + 1.5]);
title(gx, 'OneWeb Beams -3dB footprints');
drawnow;
pause(3); % wait for zoomed tiles
fig3_filename_zoom = fullfile(fig_dir, 'oneweb_projection_zoom');
exportgraphics(fig3, [fig3_filename_zoom, '.png'], 'Resolution', 300);

%% Plot 4: 3D Surf Plot in Steering Angle Space
% Convert u, v to steering angles (degrees)
Theta_U = asind(U);
Theta_V = asind(V);
theta_u_lims = asind(u_lims);
theta_v_lims = asind(v_lims);
theta_u_center = asind(BeamGrid.u_center);
theta_v_center = asind(BeamGrid.v_center);

fig4 = figure('Name', 'OneWeb Normalized Gain', 'Position', [100, 100, 600, 400]);
surf(Theta_U, Theta_V, Max_Normalized_Map, 'EdgeColor', 'none');
colormap(gca, 'turbo');
caxis([-15, 0]);
colorbar;
hold on;

scatter3(theta_u_center, theta_v_center, zeros(1, num_beams), 50, 'k', 'filled');

xlabel('Cross-Track Angle (deg)');
ylabel('Along-Track Angle (deg)');
zlabel('Normalized Gain (dB)');
title('3D Map of Maximum Normalized Gain (Steering Angle)');
axis tight;
xlim(theta_u_lims);
ylim(theta_v_lims);
zlim([-15, 0]);
view(-30, 45);

fig4_filename = fullfile(fig_dir, 'oneweb_surf_gain_steering');
exportgraphics(fig4, [fig4_filename, '.png'], 'Resolution', 300);

%% Plot 5: 2D Contours Only in Steering Angle Space
fig5 = figure('Name', 'OneWeb Beams -3dB Contours (Steering Angle)', 'Position', [100, 100, 600, 400]);
hold on;
for b = 1:num_beams
    contour(Theta_U, Theta_V, Normalized_Gain_3D(:,:,b), [-3, -3], 'LineWidth', 1.5, 'LineColor', 'b');
end
plot(theta_u_center, theta_v_center, 'k.', 'MarkerSize', 12);

xlabel('Cross-Track Angle (deg)');
ylabel('Along-Track Angle (deg)');
title('OneWeb Beams -3dB Contours');
axis equal; axis tight;
xlim(theta_u_lims);
ylim(theta_v_lims);
grid on;

fig5_filename = fullfile(fig_dir, 'oneweb_contours_only_steering');
exportgraphics(fig5, [fig5_filename, '.png'], 'Resolution', 300);

fprintf("Validation plots generated and saved successfully to %s.\n", fig_dir);