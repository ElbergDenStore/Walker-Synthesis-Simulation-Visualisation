clear all; close all; clc;

Cfg = get_cfg(1000);


calc_link = false;
plot_results = false;
use_parallel = true;



% Nordjylland
latlim = [56.5, 58.0];
lonlim = [8.0, 11.0];

people_per_ue = 300;
Cfg.Accept_Flat_UE_array = true; Cfg.Equal_UE_area = false; % Enable specific UE positions
[Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, Total_Pop] = generate_population_based_UEs(latlim, lonlim, people_per_ue);


metrics = coverage_simulator_function(Cfg, plot_results,use_parallel,calc_link);
show_interactive = true;
save_fig = false;

% show_constellation(Cfg, show_interactive, save_fig) % might be stupid
fprintf("Worst Coverage percentage " + metrics.worst_coverage_percent);

%% 1. Pick a Time Step (e.g., when the satellite is over the center)
nT = length(metrics.SimData{1}.Time);
NumUEs = length(metrics.SimData);
TotalSats = Cfg.Total_sats;
time_vec = metrics.SimData{1}.Time;

t_idx = 1; 
Re = 6371e3;
h = Cfg.Orbit_height; % From your setup (1000km)

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

figure;
histogram(all_eta);
figure;
histogram(all_az);
figure;
scatter(all_az,all_eta)

%% 2. Convert to Projection (u, v space)
% This matches your beam generation logic
u_ues = sind(all_eta) .* cosd(all_az);
v_ues = sind(all_eta) .* sind(all_az);

%% 3. Plotting
figure('Color', 'w', 'Name', 'Satellite Nadir View');
hold on; axis equal; grid on;

%% Beam data
%% 2. Math: Calculate the Boundary
Re = 6371;              % Earth Radius (km)
h_km = Cfg.Orbit_height / 1e3;       % Cfg.Orbit_heightitude (km)
Min_Elev_deg = 20;      % Boundary limit

f = 20e9;
c= 3e8;
lambda = c/f;
G = 40;

Beamwidth_deg = sqrt(32400./(10.^(G/10))) % Approximation planar array Belanis

% Antenna Size 
% G = eta * (pi * D / lambda)^2  --> D = (lambda/pi) * sqrt(G/eta)
G_linear = 10^(G*0.1);
D_meters = (lambda / pi) * sqrt(G_linear); % can be used to state
D_cm = D_meters * 100;
D = sqrt((D_cm/2)^2*pi)
% Beamwidth_deg = 4;      % Width of each individual beam

% Calculate the Max Nadir Angle to hit exactly the Min Elevation Ring
% Derived from the Law of Sines
eta_max = asind((Re / (Re + h_km)) * cosd(Min_Elev_deg));

%% Generate Hex Grid in Direction Cosine (u, v) Space
d = Beamwidth_deg; 

% In u,v space, the spacing between beams is sin(beamwidth)
du = sind(d) / 2;
rings = ceil(sind(eta_max) / du) + 5; %plus 5 to check whats limiting it

beam_data = []; %[pitch, roll, major_axis, minor_axis, plot_width, eta]

for q = -rings:rings
    for r = -rings:rings
        u = du * sqrt(3) * (q + r/2);
        v = du * 1.5 * r;
        
        if (u^2 + v^2) <= sind(eta_max)^2
            
            % Calculate Angle from Nadir
            eta = asind(sqrt(u^2 + v^2));
            
            % Calculate elliptical dimensions
            bw_major_axis = d / (cosd(eta)); % Stretched steering direction
            bw_minor_axis = d;               % Un-stretched orthogonal direction

            % Choose which one to draw. 
            plot_width = bw_major_axis; 
            % plot_width = bw_minor_axis; 
            
            
            
            roll  = asind(v);
            pitch = asind(u / cosd(roll)); %quirk from matlab sequence rotation
            % pitch = asind(u);

            beam_data = [beam_data; pitch, roll, bw_major_axis, bw_minor_axis, plot_width, eta]; 
        end
    end
end

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