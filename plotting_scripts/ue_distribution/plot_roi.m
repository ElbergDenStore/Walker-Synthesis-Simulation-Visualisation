% PLOT_ROI  Plot the simulation region of interest (ROI) without any data overlay.
%   Renders the land map using a Lambert projection for the ROI defined in
%   run_constellation_simulation_example.m (the smaller, Europe-centred longitude range).
%
%   Edit SAVE_PATH below ('' = figures/ next to this script).

save_path = '';   % output folder for ROI_Map.png

script_dir = fileparts(mfilename('fullpath'));
if isempty(save_path)
    save_path = fullfile(script_dir, 'figures');
end

%% --- ROI definition (from run_constellation_simulation_example.m, commented-out lon range) ---
lat_lim = [54 + (35/60),  83 + (40/60)];   % [54.583°, 83.667°] N
lon_lim = [-(73 + (10/60)), 33 + (30/60)]; % [-73.167°, 33.5°]

%% --- Load land areas ---
land = shaperead('landareas.shp', 'UseGeoCoords', true);

%% --- Plot ---
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 600 500]);
axesm('lambertstd', ...
    'MapLatLimit', lat_lim, ...
    'MapLonLimit', lon_lim, ...
    'Frame', 'on', ...
    'Grid', 'on', ...
    'MeridianLabel', 'on', ...
    'ParallelLabel', 'on');
axis off;

geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', ...
    'FaceColor', [0.75 0.75 0.75], 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 0.5);

% title(sprintf('Region of Interest  (%.1f°–%.1f°N,  %.1f°–%.1f°)', ...
    % lat_lim(1), lat_lim(2), lon_lim(1), lon_lim(2)), 'FontSize', 13);

set(gca, 'FontSize', 14);

out_file = fullfile(save_path, 'ROI_Map.png');
exportgraphics(f, out_file, 'Resolution', 300);
close(f);
fprintf('ROI map saved to: %s\n', out_file);
