clearvars; close all; clc;

%% 1. Configuration & Region Selection
% Toggle between 'Nordjylland', 'Denmark', or 'Full'
REGION = 'Full3000'; 

switch REGION
    case 'Nordjylland'
        latlim = [56.5, 58.0];
        lonlim = [8.0, 11.0];
        people_per_ue = 300;
        dataFile = 'Nordjylland300.mat';
    case 'Denmark'
        latlim = [54.5, 58.0];
        lonlim = [8.0, 15.5];
        people_per_ue = 3000;
        dataFile = 'Denmark3000.mat';
    case 'Full3000'
        latlim = [54.5, 83.9]; 
        lonlim = [-60, 30.0];
        people_per_ue = 3000;
        dataFile = 'Full3000.mat';
    case 'Full30000'
        latlim = [54.5, 83.9]; 
        lonlim = [-60, 30.0];
        people_per_ue = 30000;
        dataFile = 'Full30000.mat';
end

FORCE_RERUN = false; 

if exist(dataFile, 'file') && ~FORCE_RERUN 
    fprintf('Loading %s data from %s...\n', REGION, dataFile);
    load(dataFile);
else
    fprintf('Running new simulation for %s...\n', REGION);
    Cfg = get_cfg(1000); % Assuming 1000km altitude
    calc_link = false; plot_results = false; use_parallel = true;
    Cfg.Accept_Flat_UE_array = true; Cfg.Equal_UE_area = false; 
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, Total_Pop] = generate_population_based_UEs(latlim, lonlim, people_per_ue);
    metrics = coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link);
    save(dataFile);
end

%%%%% VISUALISE UES
% --- VISUAL VERIFICATION PLOT ---
fprintf('Plotting UE distribution map...\n');

% Create a clean, white-background figure window
figure('Name',  'UE Distribution', 'Color', 'w', 'Position', [100, 100, 800, 800]);

% Create geographic axes
gx = geoaxes;

% Scatter the UEs: Size 30, Red, Filled, with a slight black edge for visibility
geoscatter(gx, Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, 30, 'r', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Lock the map view specifically to Greenland's coordinates
geolimits(gx, latlim, lonlim);

% Add a basemap to show the terrain and ice sheet 
% (Other good options: 'topographic', 'streets-light', or 'bluegreen')
geobasemap(gx, 'satellite'); 

% Add a dynamic title
Cfg.NumUEs = length(Cfg.Flat_UE_array.Lats);
title(sprintf(' UE Distribution\nTotal UEs: %d (1 per %d people)', ...
    Cfg.NumUEs, round(Total_Pop/Cfg.NumUEs)), 'FontWeight', 'bold', 'FontSize', 14);

pause(3*(Cfg.NumUEs/3000)) %to load UE distribution, 3000UEs take 3 seconds
exportgraphics(gcf, 'screenshots/UEDistribution.png', 'Resolution', 300);



%% 2. System-Wide Utilization Analysis
nT = length(metrics.SimData{1}.Time);
TotalSats = Cfg.Total_sats;
time_vec = metrics.SimData{1}.Time;
time_mins = minutes(time_vec - time_vec(1));

% Map which Sat each UE is using at each second
all_SatIDs = nan(Cfg.NumUEs, nT);
for i = 1:Cfg.NumUEs
    all_SatIDs(i, :) = metrics.SimData{i}.SatID;
end

% Calculate SatUtilization [TotalSats x nT]
SatUtilization = zeros(TotalSats, nT);
for t = 1:nT
    active_sats = all_SatIDs(:, t);
    active_sats = active_sats(~isnan(active_sats));
    if ~isempty(active_sats)
        SatUtilization(:, t) = histcounts(active_sats, 0.5:(TotalSats+0.5));
    end
end

%% 3D Surface: Constellation Utilization vs Time
figure('Color', 'w', 'Name', 'Constellation Load');
[TimeGrid, SatGrid] = meshgrid(time_mins, 1:TotalSats);

% We use 'surf' to create the 3D topology
s = surf(TimeGrid, SatGrid, SatUtilization);

% --- Styling ---
s.EdgeColor = 'none'; % Smooth look
s.FaceColor = 'interp'; 
colormap(turbo); % High contrast for load levels
view(-35, 60);   % Isometric perspective
grid on;

% Add a colorbar to show UE count
cb = colorbar;
ylabel(cb, 'Connected UEs');

xlabel('Time (minutes from start)');
ylabel('Satellite ID');
zlabel('UE Count per Satellite');
title(sprintf('%s: System-Wide Satellite Utilization', REGION));

% Export the 3D view
exportgraphics(gcf, 'screenshots/constellation_surf.png', 'Resolution', 300);

% % Find satellites that had at least 1 UE at some point
% active_sat_indices = find(max(SatUtilization, [], 2) > 0);
% % Then plot only those rows:
% imagesc(time_mins, active_sat_indices, SatUtilization(active_sat_indices, :));
% exportgraphics(gcf, 'screenshots/Least1UEconstellation_surf.png', 'Resolution', 300);
%% Plot: Constellation Utilization (Active Satellites Only)
figure('Color', 'w', 'Name', 'Active Satellites Load');

% Find satellites that had at least 1 UE at some point
active_sat_indices = find(max(SatUtilization, [], 2) > 0);

% Plot only the active rows
% Note: using 1:length(active_sat_indices) on the Y-axis can be clearer 
% if the IDs are non-sequential, but usually we plot against the ID.
imagesc(time_mins, active_sat_indices, SatUtilization(active_sat_indices, :));

% --- Axis Titles & Styling ---
xlabel('Time (minutes from start)');
ylabel('Satellite ID (Constellation Index)');
title(sprintf('%s: Utilization of Active Satellites', REGION));

% Add colorbar and label it
cb = colorbar;
ylabel(cb, 'Number of Connected UEs');

% Force the Y-axis to show actual IDs and not just a scaled range
set(gca, 'YDir', 'normal'); 
grid on;

% Export
exportgraphics(gcf, 'screenshots/ActiveSatUtilization.png', 'Resolution', 300);




% --- Find the Peak Satellite and Time ---
[max_val, linear_idx] = max(SatUtilization(:));
[peakSat, t_peak] = ind2sub(size(SatUtilization), linear_idx);

fprintf('Peak Found: Sat %d at %s with %d UEs.\n', peakSat, datestr(time_vec(t_peak)), max_val);

%% 3. Geometry & Beam Grid
Re = 6371; h = Cfg.Orbit_height/1000; Min_Elev_deg = 20;
f = 20e9; c = 3e8; lambda = c/f; G = 40;
Beamwidth_deg = sqrt(32400./(10.^(G/10)));
r_beam = sind(Beamwidth_deg / 2);
eta_max = asind((Re / (Re + h)) * cosd(Min_Elev_deg));
du = sind(Beamwidth_deg) / 2;
rings = ceil(sind(eta_max) / du) + 2;
% Add arrays to store the hex grid logical coordinates
b_u = []; b_v = []; 
b_q = []; b_r = []; % NEW: To store grid logic

for q = -rings:rings
    for r = -rings:rings
        u_val = du * sqrt(3) * (q + r/2); 
        v_val = du * 1.5 * r;
        if (u_val^2 + v_val^2) <= sind(eta_max)^2
            b_u = [b_u; u_val]; 
            b_v = [b_v; v_val];
            b_q = [b_q; q]; % Save q-index
            b_r = [b_r; r]; % Save r-index
        end
    end
end

% --- THE MAGIC FIX FOR THE WATERFALL ---
% Sort the beams geographically (Row by Row, Top to Bottom, Left to Right)
% sortrows sorts primarily by the first column (v), then by the second (u)
[~, sort_idx] = sortrows([b_v, b_u], 'descend');

b_u = b_u(sort_idx);
b_v = b_v(sort_idx);
b_q = b_q(sort_idx);
b_r = b_r(sort_idx);

num_beams = length(b_u);
%% 4. Vectorize UE Projections (FILTERED BY PEAK SAT)
% We only care about UEs when they are specifically talking to peakSat
u_ues_peak = nan(Cfg.NumUEs, nT);
v_ues_peak = nan(Cfg.NumUEs, nT);

for i = 1:Cfg.NumUEs
    % Only take time steps where this UE is connected to our chosen satellite
    is_connected = (all_SatIDs(i, :) == peakSat);
    
    if any(is_connected)
        el = metrics.SimData{i}.Elevation_deg(is_connected);
        az = metrics.SimData{i}.Azimuth_deg(is_connected);
        eta_vec = asind((Re / (Re + h)) * cosd(el));
        
        u_ues_peak(i, is_connected) = sind(eta_vec) .* cosd(az + 180);
        v_ues_peak(i, is_connected) = sind(eta_vec) .* sind(az + 180);
    end
end

%% 5. Resource Analysis for the Single Peak Satellite
beam_activity_matrix = false(num_beams, nT);
for t = 1:nT
    u_t = u_ues_peak(~isnan(u_ues_peak(:,t)), t);
    v_t = v_ues_peak(~isnan(v_ues_peak(:,t)), t);
    if ~isempty(u_t)
        for b = 1:num_beams
            if any((u_t - b_u(b)).^2 + (v_t - b_v(b)).^2 <= r_beam^2)
                beam_activity_matrix(b, t) = true;
            end
        end
    end
end

%% 6. Snapshot at PEAK UTILIZATION (in Degrees)
fig1 = figure('Color', 'w', 'Visible', 'off'); hold on; axis equal; grid on;
theta = linspace(0, 2*pi, 100);

% Plot beams in degrees
b_eta = asind(sqrt(b_u.^2 + b_v.^2)); b_az = atan2d(b_v, b_u);
for b = 1:num_beams
    plot(b_eta(b)*cosd(b_az(b)) + (Beamwidth_deg/2)*cos(theta), ...
         b_eta(b)*sind(b_az(b)) + (Beamwidth_deg/2)*sin(theta), 'Color', [0.85 0.85 0.85]);
end

% Plot UEs connected to peakSat at t_peak
u_snap = u_ues_peak(:, t_peak); v_snap = v_ues_peak(:, t_peak);
ue_eta = asind(sqrt(u_snap.^2 + v_snap.^2)); ue_az = atan2d(v_snap, u_snap);
scatter(ue_eta.*cosd(ue_az), ue_eta.*sind(ue_az), 3, 'r');



plot(eta_max*cos(theta), eta_max*sin(theta), 'k--', 'LineWidth', 2);
title(sprintf('Sat %d Peak Load (%d/%d UEs)\nTime: %s', peakSat, max_val, Cfg.NumUEs, datestr(time_vec(t_peak))));
xlabel('Degrees off Nadir'); ylabel('Degrees off Nadir');
exportgraphics(fig1, 'screenshots/peak_sat_fov.png', 'Resolution', 300);



% % --- Plot Waterfall for this satellite only ---
% fig2 = figure('Color', 'w', 'Visible', 'off');
% imagesc(time_mins, 1:num_beams, beam_activity_matrix);
% colormap([1 1 1; 0.6 0.2 0.2]); % Redish for this specific sat's load
% xlabel('Time (mins)'); ylabel('Beam ID');
% title(sprintf('Waterfall: Switching Profile for Sat %d', peakSat));
% exportgraphics(fig2, 'screenshots/peak_sat_waterfall.png', 'Resolution', 300);


% fig2b = figure('Color', 'w', 'Visible', 'off');
% imagesc(time_mins, 1:num_beams, beam_activity_matrix);
% colormap([1 1 1; 0.6 0.2 0.2]); % Redish for this specific sat's load
% xlabel('Time (mins)'); ylabel('Beam ID');
% ylim([500 600])
% xlim([time_mins(1) time_mins(round(nT/10))]);
% title(sprintf('Zoomed Waterfall: Switching Profile for Sat %d', peakSat));
% exportgraphics(fig2b, 'screenshots/ZoomPeak_sat_waterfall.png', 'Resolution', 300);

fprintf('Done! Saved peak analysis for Sat %d.\n', peakSat);


%% Resource Load Analysis & Beam Assignment (Single Sat Perspective)
active_beams_time = zeros(nT, 1);
max_ues_per_beam_time = zeros(nT, 1);
beam_ue_count_matrix = zeros(num_beams, nT); % Preallocate here!
r_beam_sq = r_beam^2;

fprintf('Calculating resource load for Peak Sat %d...\n', peakSat);

for t = 1:nT
    u_t = u_ues_peak(~isnan(u_ues_peak(:,t)), t);
    v_t = v_ues_peak(~isnan(v_ues_peak(:,t)), t);
    
    if ~isempty(u_t)
        % 1. Calculate distance from every active UE to every beam center
        dist_sq_matrix = (u_t - b_u').^2 + (v_t - b_v').^2;
        
        % 2. Find the index of the closest beam for each UE (No double counting!)
        [min_dist_sq, best_beam_idx] = min(dist_sq_matrix, [], 2);
        
        % 3. Filter: Only count if the UE is actually within the beam radius
        % valid_connection = min_dist_sq <= r_beam_sq;
        % best_beam_idx = best_beam_idx(valid_connection);
        
        if ~isempty(best_beam_idx)
            % 4. Count UEs per beam for this time step
            temp_beam_counts = histcounts(best_beam_idx, 0.5:(num_beams+0.5));
            
            % 5. Metrics
            active_beams_time(t) = sum(temp_beam_counts > 0);
            max_ues_per_beam_time(t) = max(temp_beam_counts);
            
            % Save to matrix for Waterfall and Heatmaps
            beam_ue_count_matrix(:, t) = temp_beam_counts';
        end
    end
end

%% 2. Plotting: The Resource Load (Active Beams)
figA = figure('Color', 'w', 'Position', [100 100 900 450]);
tiledlayout(2,1, 'TileSpacing', 'compact');

% Top Plot: Total Beams Needed
nexttile;
plot(time_mins, active_beams_time, 'LineWidth', 2, 'Color', [0 0.447 0.741]);
grid on; ylabel('Active Beams');
title(sprintf('Sat %d Resource Profile: Active Beams', peakSat));
yline(mean(active_beams_time(active_beams_time>0)), '--k', 'Avg Active');

% Bottom Plot: Max Congestion (UEs/Beam)
nexttile;
plot(time_mins, max_ues_per_beam_time, 'LineWidth', 2, 'Color', [0.85 0.32 0.1]);
grid on; ylabel('Max UEs in 1 Beam');
xlabel('Time (minutes from start)');
title('Beam Congestion Level');

exportgraphics(figA, 'screenshots/peak_sat_resource_load.png', 'Resolution', 300);

%% 3D Waterfall: Spatial Congestion vs Time
figure('Color', 'w', 'Name', 'Beam Load Waterfall');

% Plotting as a 3D surface or Heatmap (Using the correctly populated matrix)
imagesc(time_mins, 1:num_beams, beam_ue_count_matrix);
colormap(parula); % Blue to Yellow
cb = colorbar;
ylabel(cb, 'UEs per Beam');
xlabel('Time (mins)'); ylabel('Beam ID');
title(sprintf('Sat %d: Spatial Congestion Waterfall', peakSat));

exportgraphics(gcf, 'screenshots/peak_sat_congestion_waterfall.png', 'Resolution', 300);

%% 7. Spatial Load Heatmap (Steering Angle Perspective)

% 1. Extract the UE counts for all beams at the peak time
counts_at_peak = beam_ue_count_matrix(:, t_peak);
active_beam_indices = find(counts_at_peak > 0);
active_counts = counts_at_peak(active_beam_indices);

% 2. Calculate the 90th Percentile for the color scale
if isempty(active_counts)
    cap_val = 1; % Fallback
else
    cap_val = prctile(active_counts, 90);
    if cap_val == 0; cap_val = max(active_counts); end % Prevent 0-cap 
end
fprintf('Capping color scale at %d UEs (90th Percentile).\n', round(cap_val));

% --- CONVERT TO DEGREES ---
% Map u,v to Nadir Angle (eta) and Azimuth (az)
b_eta = asind(sqrt(b_u.^2 + b_v.^2));
b_az  = atan2d(b_v, b_u);

% Convert to Cartesian X/Y where distance = degrees off nadir
b_x_deg = b_eta .* cosd(b_az);
b_y_deg = b_eta .* sind(b_az);
r_beam_deg = Beamwidth_deg / 2; % Beam radius is now simply half the HPBW in degrees

% 3. Prepare coordinates for drawing ELLIPSES in Degree-Space
theta = linspace(0, 2*pi, 40); 
X_all = zeros(length(theta), num_beams);
Y_all = zeros(length(theta), num_beams);

for b = 1:num_beams
    eta = b_eta(b);
    az  = b_az(b);
    
    % Prevent division by zero for the exact nadir beam
    if eta < 1e-3
        map_distortion_factor = 1;
    else
        % Convert eta to radians for the numerator, use sind for the denominator
        map_distortion_factor = (eta * pi/180) / sind(eta);
    end
    
    % The major axis stretches based on physical array Scan Loss
    r_major = r_beam_deg / cosd(eta); 
    
    % The minor axis stretches based on the 2D plot projection distortion
    r_minor = r_beam_deg * map_distortion_factor; 
    
    % Unrotated ellipse centered at 0,0
    x_local = r_major * cos(theta);
    y_local = r_minor * sin(theta);
    
    % Rotate the ellipse so the stretch aligns with the steering azimuth
    x_rot = x_local * cosd(az) - y_local * sind(az);
    y_rot = x_local * sind(az) + y_local * cosd(az);
    
    % Translate to the beam center
    X_all(:, b) = b_x_deg(b) + x_rot;
    Y_all(:, b) = b_y_deg(b) + y_rot;
end

% Easily grab the active ones using the indices
X_active = X_all(:, active_beam_indices);
Y_active = Y_all(:, active_beam_indices);

% =========================================================================
% PLOT A: FULL FOV SPATIAL LOAD (IN DEGREES)
% =========================================================================
fig_spatial = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 800 800]);
hold on; axis equal; box on;

% Draw the inactive beams as faint outlines
plot(X_all, Y_all, 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5);

% Draw the active beams as filled patches colored by UE count
patch(X_active, Y_active, active_counts', 'EdgeColor', 'k', 'LineWidth', 0.5, 'FaceAlpha', 0.85);

% Plot the Horizon Ring (in degrees)
plot(eta_max*cos(theta), eta_max*sin(theta), 'k--', 'LineWidth', 2);

% Styling and Colorbar
colormap(turbo); 
try clim([0 cap_val]); catch; caxis([0 cap_val]); end % Compatibility for older MATLAB
cb = colorbar;
ylabel(cb, sprintf('Connected UEs (Capped at 90th pct: %d)', round(cap_val)));

xlabel('X Steering Angle (Degrees off Nadir)'); 
ylabel('Y Steering Angle (Degrees off Nadir)');
title(sprintf('Sat %d: FoV Beam Load Heatmap\nTime: %s', peakSat, datestr(time_vec(t_peak))));
xlim([-eta_max-2, eta_max+2]); ylim([-eta_max-2, eta_max+2]);

exportgraphics(fig_spatial, 'screenshots/peak_sat_spatial_full_deg.png', 'Resolution', 300);

% =========================================================================
% PLOT B: ZOOMED IN "FREQUENCY REUSE" VIEW (IN DEGREES)
% =========================================================================

% 1. Find the index of the beam with the absolute maximum UE count
[max_ue_val, max_idx] = max(active_counts);

% 2. Map that index back to the original Beam ID
b_idx_max = active_beam_indices(max_idx);

% 3. Extract the exact center coordinates of that specific beam
center_x_deg = b_x_deg(b_idx_max);
center_y_deg = b_y_deg(b_idx_max);

fprintf('Zooming in on the busiest beam: ID %d with %d UEs.\n', b_idx_max, max_ue_val);

fig_zoom = figure('Color', 'w', 'Visible', 'off', 'Position', [150 150 900 900]);
hold on; axis equal; box on; grid on;

% Draw all beams (gray outlines)
plot(X_all, Y_all, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);

% Draw active beams (filled)
patch(X_active, Y_active, active_counts', 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.85);

% --- ADD LABELS FOR REUSE PLANNING ---
for i = 1:length(active_beam_indices)
    b_idx = active_beam_indices(i);
    text(b_x_deg(b_idx), b_y_deg(b_idx), sprintf('ID: %d\nUEs: %d', b_idx, counts_at_peak(b_idx)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'w');
end

% --- Highlight the absolute max beam with a thicker outline ---
eta_max_beam = b_eta(b_idx_max);
az_max_beam  = b_az(b_idx_max);

% 1. Calculate the map distortion factor for the peak beam
if eta_max_beam < 1e-3
    map_dist_max = 1;
else
    map_dist_max = (eta_max_beam * pi/180) / sind(eta_max_beam);
end

% 2. Apply both Scan Loss (Major) and Projection Distortion (Minor)
r_major_max = r_beam_deg / cosd(eta_max_beam);
r_minor_max = r_beam_deg * map_dist_max;

% 3. Build local ellipse
x_loc_max = r_major_max * cos(theta);
y_loc_max = r_minor_max * sin(theta);

% 4. Rotate to match steering azimuth
x_rot_max = x_loc_max * cosd(az_max_beam) - y_loc_max * sind(az_max_beam);
y_rot_max = x_loc_max * sind(az_max_beam) + y_loc_max * cosd(az_max_beam);

% 5. Draw the red ring exactly over the patch
plot(center_x_deg + x_rot_max, center_y_deg + y_rot_max, 'r', 'LineWidth', 3);

colormap(turbo); 
try clim([0 cap_val]); catch; caxis([0 cap_val]); end
cb_zoom = colorbar; ylabel(cb_zoom, 'Connected UEs');

% 4. Apply the zoom limits centered directly on the busiest beam
zoom_radius = 5 * (r_beam_deg * 2); 
xlim([center_x_deg - zoom_radius, center_x_deg + zoom_radius]);
ylim([center_y_deg - zoom_radius, center_y_deg + zoom_radius]);

xlabel('X Steering Angle (deg)'); ylabel('Y Steering Angle (deg)');
title(sprintf('Zoomed Active Region: Sat %d\nCentered on Peak Beam %d (Max Load)', peakSat, b_idx_max));

exportgraphics(fig_zoom, 'screenshots/peak_sat_spatial_zoomed_deg.png', 'Resolution', 300);

%% 8. Beam Load Percentile Distribution

% 1. Extract and sort the active beam counts
counts_at_peak = beam_ue_count_matrix(:, t_peak);
active_counts = counts_at_peak(counts_at_peak > 0);

% Sort from least loaded to most loaded
sorted_counts = sort(active_counts, 'ascend');

% 2. Create the Percentile X-axis (0% to 100%)
num_active = length(sorted_counts);
if num_active > 1
    percentiles = linspace(0, 100, num_active);
else
    percentiles = 100; % Edge case if only 1 beam is active
end

% 3. Calculate total connected UEs for the title
total_ues_connected = sum(sorted_counts);

% 4. Plotting
fig_hist = figure('Color', 'w', 'Visible', 'off', 'Position', [200 200 800 500]);

% Using 'area' for a clean, filled histogram look
area(percentiles, sorted_counts, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', [0 0.3 0.6], 'LineWidth', 1.5);
hold on; grid on;

% Add helpful statistical reference lines
mean_load = mean(sorted_counts);
p90_load = prctile(sorted_counts, 90);

yline(mean_load, '--r', sprintf('Mean Load (%.1f UEs)', mean_load), 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
yline(p90_load, '--k', sprintf('90th Percentile (%.1f UEs)', p90_load), 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');

% Formatting
xlabel('Active Beam Percentile (%)');
ylabel('Number of Connected UEs');
title(sprintf('Sat %d: Beam Load Distribution\nTime: %s | Total Connected UEs: %d (Across %d Beams)', ...
    peakSat, datestr(time_vec(t_peak)), total_ues_connected, num_active));

% Set limits for a clean look
xlim([0 100]);
ylim([0 max(sorted_counts) * 1.1]); % Add 10% headroom to the top

exportgraphics(fig_hist, 'screenshots/peak_sat_load_percentile.png', 'Resolution', 300);
fprintf('Percentile Load Distribution plot generated! Check screenshots folder.\n');


%% 8. Beam Load CDF (Cumulative Distribution Function)

% 1. Extract and sort the active beam counts
counts_at_peak = beam_ue_count_matrix(:, t_peak);
active_counts = counts_at_peak(counts_at_peak > 0);
sorted_counts = sort(active_counts, 'ascend');

% 2. Calculate the Empirical CDF (0 to 100%)
num_active = length(sorted_counts);
if num_active > 0
    percentiles = (1:num_active) / num_active * 100;
else
    percentiles = 0;
end

% 3. Calculate total connected UEs for the title
total_ues_connected = sum(sorted_counts);

% 4. Plotting the CDF
fig_cdf = figure('Color', 'w', 'Visible', 'off', 'Position', [200 200 800 500]);

% Plot the main CDF curve
plot(sorted_counts, percentiles, '-', 'LineWidth', 3, 'Color', [0 0.447 0.741]);
hold on; grid on;

% Optional: Add a light fill under the curve for a polished aesthetic
area(sorted_counts, percentiles, 'FaceColor', [0 0.447 0.741], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% 5. Add helpful statistical reference lines on the X-axis
mean_load = mean(sorted_counts);
p90_load = prctile(sorted_counts, 90);

xline(mean_load, '--r', sprintf('Mean (%.1f UEs)', mean_load), ...
    'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 10);
    
xline(p90_load, '--k', sprintf('90th Percentile (%.1f UEs)', p90_load), ...
    'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 10);

% 6. Formatting
xlabel('Number of Connected UEs per Beam');
ylabel('Cumulative Percentage of Beams (%)');
title(sprintf('Sat %d: Beam Load CDF\nTime: %s | Total Connected UEs: %d/%d (Across %d Beams)', ...
    peakSat, datestr(time_vec(t_peak)), total_ues_connected,Cfg.NumUEs, num_active));

% Set axis limits for a clean look
if max(sorted_counts) > 0
    xlim([0 max(sorted_counts) * 1.05]); % Add 5% headroom to the right
else
    xlim([0 1]);
end
ylim([0 100]);

exportgraphics(fig_cdf, 'screenshots/peak_sat_load_cdf.png', 'Resolution', 300);
fprintf('CDF Load Distribution plot generated! Check screenshots folder.\n');