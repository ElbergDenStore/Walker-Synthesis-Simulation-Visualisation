% 1. Define your target coordinates
target_lat = 57.05;
target_lon = 9.92;

% 2. Read the GeoTIFF file
% (Make sure to replace this with your actual downloaded filename)
filename = 'ppp_2020_1km_Aggregated.tif';
[pop_data, R] = readgeoraster(filename);

% 3. Define the geographic bounding box for Denmark
% Denmark roughly spans: Lat 54.5°N to 58.0°N, Lon 8.0°E to 15.5°E
latlim = [54.5, 58.0];
lonlim = [8.0, 15.5];

% 4. Crop the dataset down to just the Denmark region 
% This saves memory and prevents the figure from crashing
fprintf('Cropping data to Denmark...\n');
[denmark_data, denmark_R] = geocrop(pop_data, R, latlim, lonlim);

% 5. Clean up the data
% Convert to double for plotting and handle 'NoData' values (often negative in these datasets)
denmark_data = double(denmark_data);
denmark_data(denmark_data < 0) = NaN; 

% 6. Create the map figure
figure('Name', 'Population Density of Denmark', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Initialize a geographic map axis centered on our limits
worldmap(latlim, lonlim);
title('Population Density of Denmark (People per pixel)', 'FontWeight', 'bold');

% 7. Plot the cropped raster data onto the map
% Using 'texturemap' or 'surface' drapes the data over the coordinates
geoshow(denmark_data, denmark_R, 'DisplayType', 'surface');

% 8. Make the plot visually readable
% Population data is highly skewed (cities are massive spikes), so we cap the color axis
caxis([0, 500]); % Adjust this 500 up or down depending on how dense Copenhagen looks
colormap(parula); % Standard MATLAB colormap
cb = colorbar;
cb.Label.String = 'Estimated People per Pixel';

% 9. Overlay standard coastlines so you can actually see the shape of the country clearly
load coastlines;
plotm(coastlat, coastlon, 'Color', 'white', 'LineWidth', 1.5);

fprintf('Map plotted successfully!\n');



% 3. Define the geographic bounding box for Greenland
% Greenland spans from roughly 59°N up to 83.6°N, and 74°W to 11°W
latlim = [59.0, 83.9];
lonlim = [-74.0, -11.0];

% 4. Crop the dataset down to the Greenland region
fprintf('Cropping data to Greenland...\n');
[greenland_data, greenland_R] = geocrop(pop_data, R, latlim, lonlim);

% 5. Clean up the data
% Convert to double and set ocean/ice 'NoData' values to NaN so they don't plot
greenland_data = double(greenland_data);
greenland_data(greenland_data < 0) = NaN; 

% 6. Create the map figure
figure('Name', 'Population Density of Greenland', 'Color', 'w', 'Position', [100, 100, 800, 700]);

% Initialize a geographic map axis centered on Greenland
worldmap(latlim, lonlim);
title('Population Density of Greenland', 'FontWeight', 'bold');

% 7. Plot the cropped raster data onto the map
geoshow(greenland_data, greenland_R, 'DisplayType', 'surface');

% 8. Make the plot visually readable
% Greenland is highly unpopulated except for a few coastal towns (like Nuuk).
% We cap the color axis at 50 so the smaller settlements actually light up.
caxis([0, 50]); 
colormap(parula);
cb = colorbar;
cb.Label.String = 'Estimated People per Pixel';

% 9. Overlay standard coastlines to give the data shape
load coastlines;
plotm(coastlat, coastlon, 'Color', 'white', 'LineWidth', 1.5);

fprintf('Map plotted successfully!\n');


% --- VISUAL VERIFICATION PLOT ---
fprintf('Plotting UE distribution map...\n');

% Create a clean, white-background figure window
figure('Name', 'Greenland UE Distribution', 'Color', 'w', 'Position', [100, 100, 800, 800]);

% Create geographic axes
gx = geoaxes;

% Scatter the UEs: Size 30, Red, Filled, with a slight black edge for visibility
geoscatter(gx, UE_lats_flat, UE_lons_flat, 30, 'r', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Lock the map view specifically to Greenland's coordinates
geolimits(gx, [59.0, 84.0], [-74.0, -11.0]);

% Add a basemap to show the terrain and ice sheet 
% (Other good options: 'topographic', 'streets-light', or 'bluegreen')
geobasemap(gx, 'satellite'); 

% Add a dynamic title
title(sprintf('Greenland UE Distribution\nTotal UEs: %d (1 per %d people)', ...
    length(UE_lats_flat), people_per_ue), 'FontWeight', 'bold', 'FontSize', 14);

fprintf('Close the map figure to continue with the simulation, or just let it run in the background.\n');
% ---------------------------------