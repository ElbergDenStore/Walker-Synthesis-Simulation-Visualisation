clear all; close all; clc;

Cfg = get_cfg(1000);


calc_link = false;
plot_results = false;
use_parallel = true;

metrics = coverage_simulator_function(Cfg, plot_results,use_parallel,calc_link);
show_interactive = true;
save_fig = false;

show_constellation(Cfg, show_interactive, save_fig) % might be stupid
fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);

% Nordjylland
latlim = [56.5, 58.0];
lonlim = [8.0, 11.0];

people_per_ue = 300;
Cfg.Accept_Flat_UE_array = true; Cfg.Equal_UE_area = false; % Enable specific UE positions
[Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, Total_Pop] = generate_population_based_UEs(latlim, lonlim, people_per_ue);


%% 1. Pick a Time Step (e.g., when the satellite is over the center)
t_idx = round(nT / 2); 
Re = 6371e3;
h = alt; % From your setup (1000km)

% Collect data for all UEs at this timestamp
all_eta = [];
all_az  = [];

for i = 1:NumUEs
    el = metrics.SimData{i}.Elevation_deg(t_idx);
    az = metrics.SimData{i}.Azimuth_deg(t_idx);
    
    if ~isnan(el)
        % Convert Elevation to Nadir Angle
        eta = asind((Re / (Re + h)) * cosd(el));
        
        all_eta = [all_eta; eta];
        % Azimuth from Sat to UE is approx Azimuth from UE to Sat + 180
        all_az  = [all_az; mod(az + 180, 360)]; 
    end
end

%% 2. Convert to Projection (u, v space)
% This matches your beam generation logic
u_ues = sind(all_eta) .* cosd(all_az);
v_ues = sind(all_eta) .* sind(all_az);

%% 3. Plotting
figure('Color', 'w', 'Name', 'Satellite Nadir View');
hold on; axis equal; grid on;

% Plot the Beam Boundaries (Circles)
theta = linspace(0, 2*pi, 50);
for b = 1:size(beam_data, 1)
    % Simplified: Plot beam circles using their centers (u,v) and width
    b_u = sind(beam_data(b,6)) * cosd(beam_data(b,1)); % rough center
    b_v = sind(beam_data(b,6)) * sind(beam_data(b,2));
    
    % Using your Beamwidth_deg transformed to u-v scale
    r_beam = sind(Beamwidth_deg/2);
    plot(sind(beam_data(b,6))*cos(theta) + u_ues(1)*0, ... % This is a placeholder
         'Color', [0.8 0.8 0.8]); % Light gray for beams
end

% Plot the UEs
scatter(u_ues, v_ues, 15, 'filled', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5);

% Plot the Min Elevation Ring
circle_u = sind(eta_max) * cos(theta);
circle_v = sind(eta_max) * sind(theta);
plot(circle_u, circle_v, 'k--', 'LineWidth', 2, 'DisplayName', 'Coverage Edge');

xlabel('u (Direction Cosine)');
ylabel('v (Direction Cosine)');
title(sprintf('Snapshot: UE Distribution in Satellite FoV\nTime: %s', string(time_vec(t_idx))));