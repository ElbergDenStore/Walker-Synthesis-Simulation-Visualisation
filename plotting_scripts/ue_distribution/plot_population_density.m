clearvars; close all; clc;

% Population density map for the same region as the Jutland zoom used in
% plotting_scripts/beams_on_earth_matlab_viewer.m:
% geolimits([lat-0.8, lat+0.8], [lon-1.5, lon+1.5]) with lat=57, lon=9.3

center_lat = 57.0;
center_lon = 9.3;
latlim = [center_lat - 0.8, center_lat + 0.8];
lonlim = [center_lon - 1.5, center_lon + 1.5];
%%% FILE IS HUGE AND CAN BE DOWNLOADED FROM https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/0_Mosaicked/ppp_2020_1km_Aggregated.tif
filename = 'ppp_2020_1km_Aggregated.tif';
script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..');  % matlab_code/
path_setup();

candidate_files = {
    % fullfile(repo_root, filename)
    fullfile(repo_root, 'functions/data', filename)
    % fullfile(script_dir, filename)
    % fullfile(pwd, filename)
};

tif_path = '';
for i = 1:numel(candidate_files)
    if isfile(candidate_files{i})
        tif_path = candidate_files{i};
        break;
    end
end

if strlength(string(tif_path)) == 0
    error(['Could not find %s. Checked:\n - %s\n - %s\n - %s\n - %s'], filename, ...
        candidate_files{1}, candidate_files{2}, candidate_files{3}, candidate_files{4});
end

fprintf('Loading raster: %s\n', tif_path);
[pop_data, R] = readgeoraster(tif_path);

fprintf('Cropping raster to Jutland zoom region...\n');
[zoom_data, zoom_R] = geocrop(pop_data, R, latlim, lonlim);

% Convert to double and clean no-data values.
zoom_data = double(zoom_data);
zoom_data(zoom_data < 0) = NaN;

% Convert people per pixel -> people per km^2
% Pixel area varies with latitude for geographic rasters.
[n_rows, n_cols] = size(zoom_data);
[col_grid, row_grid] = meshgrid(1:n_cols, 1:n_rows);
[lat_center, ~] = intrinsicToGeographic(zoom_R, col_grid, row_grid);

dlat_deg = zoom_R.CellExtentInLatitude;
dlon_deg = zoom_R.CellExtentInLongitude;

km_per_deg_lat = 111.32;
km_per_deg_lon = 111.32 .* cosd(lat_center);
pixel_area_km2 = (km_per_deg_lat * dlat_deg) .* (km_per_deg_lon * dlon_deg);

pop_density_km2 = zoom_data ./ pixel_area_km2;
pop_density_km2(~isfinite(pop_density_km2)) = NaN;

fig = figure('Name', 'Population Density - Jutland Zoom', 'Color', 'w', 'Position', [100, 100, 600, 400]);
worldmap(latlim, lonlim);
geoshow(pop_density_km2, zoom_R, 'DisplayType', 'surface');
load coastlines;
plotm(coastlat, coastlon, 'w', 'LineWidth', 1.0);

% Cap color scale at the 99th percentile to avoid city outliers dominating.
vals = pop_density_km2(isfinite(pop_density_km2));
if isempty(vals)
    cmax = 1;
else
    cmax = prctile(vals, 99);
end
caxis([0, cmax]);
colormap(turbo);
cb = colorbar;
cb.Label.String = 'People per km^2';

title(sprintf('Population Density (people/km^2) | [%.1f to %.1f N, %.1f to %.1f E]', ...
    latlim(1), latlim(2), lonlim(1), lonlim(2)), 'FontWeight', 'bold');

out_dir = fullfile(repo_root, 'figures', 'population_density');
if ~isfolder(out_dir)
    mkdir(out_dir);
end
out_file = fullfile(out_dir, 'jutland_population_density_km2.png');
exportgraphics(fig, out_file, 'Resolution', 300);

fprintf('Saved: %s\n', out_file);