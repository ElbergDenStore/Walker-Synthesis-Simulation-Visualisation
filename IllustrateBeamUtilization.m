clearvars; close all; clc;

%% 1. Configuration & Region Selection
% Toggle between 'Nordjylland', 'Denmark', or 'Full'
REGION = 'Full'; 

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
    case 'Full'
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
pause(3)
exportgraphics(gcf, 'screenshots/UEDistribution.png', 'Resolution', 300);



%% 2. System-Wide Utilization Analysis
nT = length(metrics.SimData{1}.Time);
NumUEs = length(metrics.SimData);
TotalSats = Cfg.Total_sats;
time_vec = metrics.SimData{1}.Time;
time_mins = minutes(time_vec - time_vec(1));

% Map which Sat each UE is using at each second
all_SatIDs = nan(NumUEs, nT);
for i = 1:NumUEs
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
b_u = []; b_v = []; 
for q = -rings:rings
    for r = -rings:rings
        u_val = du * sqrt(3) * (q + r/2); v_val = du * 1.5 * r;
        if (u_val^2 + v_val^2) <= sind(eta_max)^2
            b_u = [b_u; u_val]; b_v = [b_v; v_val];
        end
    end
end
num_beams = length(b_u);

%% 4. Vectorize UE Projections (FILTERED BY PEAK SAT)
% We only care about UEs when they are specifically talking to peakSat
u_ues_peak = nan(NumUEs, nT);
v_ues_peak = nan(NumUEs, nT);

for i = 1:NumUEs
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
scatter(ue_eta.*cosd(ue_az), ue_eta.*sind(ue_az), 3, 'r', 'filled');



plot(eta_max*cos(theta), eta_max*sin(theta), 'k--', 'LineWidth', 2);
title(sprintf('Sat %d Peak Load (%d/%d UEs)\nTime: %s', peakSat, max_val, Cfg.NumUEs, datestr(time_vec(t_peak))));
xlabel('Degrees off Nadir'); ylabel('Degrees off Nadir');
exportgraphics(fig1, 'screenshots/peak_sat_fov.png', 'Resolution', 300);

% --- Plot Waterfall for this satellite only ---
fig2 = figure('Color', 'w', 'Visible', 'off');
imagesc(time_mins, 1:num_beams, beam_activity_matrix);
colormap([1 1 1; 0.6 0.2 0.2]); % Redish for this specific sat's load
xlabel('Time (mins)'); ylabel('Beam ID');
title(sprintf('Waterfall: Switching Profile for Sat %d', peakSat));
exportgraphics(fig2, 'screenshots/peak_sat_waterfall.png', 'Resolution', 300);

fprintf('Done! Saved peak analysis for Sat %d.\n', peakSat);


%% Resource Load Analysis (Single Sat Perspective)
active_beams_time = zeros(nT, 1);
max_ues_per_beam_time = zeros(nT, 1);
r_beam_sq = r_beam^2;

fprintf('Calculating resource load for Peak Sat %d...\n', peakSat);

for t = 1:nT
    % Get only UEs connected to our peakSat right now
    u_t = u_ues_peak(~isnan(u_ues_peak(:,t)), t);
    v_t = v_ues_peak(~isnan(v_ues_peak(:,t)), t);
    
    if ~isempty(u_t)
        temp_beam_counts = zeros(num_beams, 1);
        for b = 1:num_beams
            % Count how many UEs fall into THIS specific beam
            dists_sq = (u_t - b_u(b)).^2 + (v_t - b_v(b)).^2;
            temp_beam_counts(b) = sum(dists_sq <= r_beam_sq);
        end
        
        % Metrics for this second
        active_beams_time(t) = sum(temp_beam_counts > 0);
        max_ues_per_beam_time(t) = max(temp_beam_counts);
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

% Create a matrix of [Beams x Time] containing UE counts
beam_ue_count_matrix = zeros(num_beams, nT);
for t = 1:nT
    u_t = u_ues_peak(~isnan(u_ues_peak(:,t)), t);
    v_t = v_ues_peak(~isnan(v_ues_peak(:,t)), t);
    if ~isempty(u_t)
        for b = 1:num_beams
            beam_ue_count_matrix(b,t) = sum((u_t - b_u(b)).^2 + (v_t - b_v(b)).^2 <= r_beam_sq);
        end
    end
end

% Plotting as a 3D surface or Heatmap
imagesc(time_mins, 1:num_beams, beam_ue_count_matrix);
colormap(parula); % Blue to Yellow
cb = colorbar;
ylabel(cb, 'UEs per Beam');
xlabel('Time (mins)'); ylabel('Beam ID');
title(sprintf('Sat %d: Spatial Congestion Waterfall', peakSat));

exportgraphics(gcf, 'screenshots/peak_sat_congestion_waterfall.png', 'Resolution', 300);