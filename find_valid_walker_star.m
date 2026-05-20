clear; close all; clc;

heights_km            = 500:0.1:1200;
Lat_range_deg         = [54+(35/60), 83+(40/60)]; %54°35N Denmark minimum, 83°40N Greenland max
Min_elevation_UE      = 20;
minimum_lat_deg = min(Lat_range_deg);

num_steps = length(heights_km);

% Preallocate arrays to store results
total_sats= zeros(1, num_steps); planes = zeros(1, num_steps); sats = zeros(1, num_steps);

for i = 1:length(heights_km) 
    [Num_planes, Sats_per_plane, Total_sats] = calculate_walker_star(heights_km(i), minimum_lat_deg, Min_elevation_UE);
    total_sats(i)= Total_sats; planes(i) = Num_planes; sats(i) = Sats_per_plane;
end

% --- Plotting the Results ---
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
folder_name = sprintf('Walker-Star_%s', date_str);
out_dir = fullfile('simulation_output', folder_name);
            
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end


% --- Plotting the Results ---
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
folder_name = sprintf('Walker-Star_%s', date_str);
out_dir = fullfile('simulation_output', folder_name);
            
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Create figure with a white background
f1 = figure('Color', 'w', 'Position', [100, 100, 700, 450]);

% --- Left Axis (Planes & Sats/Plane) ---
yyaxis left;
% Using MATLAB's modern, colorblind-friendly hex palette with mixed line styles
plot(heights_km, planes, '-', 'Color', '#D95319', 'LineWidth', 1.5); hold on; 
plot(heights_km, sats, '-', 'Color', '#0072BD', 'LineWidth', 1.5);

% Make the left axis and label neutral (black) since it tracks two different variables
ylabel('Planes | Sats/Plane', 'Color', 'k'); 
ax = gca;
ax.YAxis(1).Color = 'k'; 
ylim([0, max(sats)]);

% --- Right Axis (Total Satellites) ---
yyaxis right;
plot(heights_km, total_sats, '-', 'Color', 'k', 'LineWidth', 2);
ylabel('Total Satellites', 'Color', 'k');
ax.YAxis(2).Color = 'k';
ylim([0, max(total_sats)*1.5]);
% --- Formatting & Aesthetics ---
grid on;
ax.GridAlpha = 0.25; % Softer grid lines to keep focus on the data
xlabel('Orbital altitude (km)');
title(sprintf('Walker-Star Configuration | Lat: %0.1f^{\\circ} | \\epsilon_{min}: %0.1f^{\\circ}', Min_elevation_UE, minimum_lat_deg));
% Legend improvements
lgd = legend('Planes', 'Sats/Plane', 'Total Satellites', 'Location', 'best');
% lgd.Box = 'off'; % Removes the distracting border around the legend

% IEEE Standard Fonts
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

% Export
exportgraphics(f1, fullfile(out_dir, 'Walker-Star.png'), 'Resolution', 300);
close(f1);


