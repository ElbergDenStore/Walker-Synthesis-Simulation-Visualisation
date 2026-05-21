function [UE_lats_flat, UE_lons_flat, Total_Pop] = generate_population_based_UEs(lat_vec, lon_vec, people_per_ue)
    % 1. Load the population data
    
    filename = 'ppp_2020_1km_Aggregated.tif';
    fprintf('Loading GeoTIFF for UE generation...\n');

    % Simply rely on the MATLAB path. If you added the functions folder using path_setup()
    % or addpath(), MATLAB will find it automatically.
    tif_path = which(filename);
    
    if isempty(tif_path)
        error('Could not find %s anywhere on the MATLAB path. Ensure the parent directory is added (e.g. by running path_setup). FILE IS HUGE AND CAN BE DOWNLOADED FROM https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/0_Mosaicked/ppp_2020_1km_Aggregated.tif.', filename);
    end

    [pop_data, R] = readgeoraster(tif_path);
    
    % 2. Get the bounding box from your simulation config
    latlim = [min(lat_vec), max(lat_vec)];
    lonlim = [min(lon_vec), max(lon_vec)];
    
    % 3. Crop the massive dataset down to your simulation area
    [local_data, local_R] = geocrop(pop_data, R, latlim, lonlim);
    
    % 4. Clean the data (Remove negative NoData values/oceans)
    local_data = double(local_data);
    local_data(local_data < 0 | isnan(local_data)) = 0;

    
    % ---------------------------------------------------------
    % 5. THE EXCLUSION MASK: Remove Canada
    % We create a 2D grid of the latitudes and longitudes for every pixel
    [cols, rows] = meshgrid(1:size(local_data,2), 1:size(local_data,1));
    [lats_grid, lons_grid] = intrinsicToGeographic(local_R, cols, rows);
    
    % Exclude Iceland (Everything south of 67.5N and east of 25W)
    % iceland_mask = (lats_grid < 67.5) & (lons_grid > -25.0);
    % local_data(iceland_mask) = 0;
    
    % Exclude Canada (Using a rough staircase cut along the western border)
    canada_mask1 = (lats_grid < 75.0) & (lons_grid < -55.0);
    canada_mask2 = (lats_grid >= 75.0) & (lons_grid < -75.0);
    local_data(canada_mask1 | canada_mask2) = 0;
    % ---------------------------------------------------------

    
    % 5. Calculate Total Population and Target Number of UEs
    Total_Pop = sum(local_data, 'all');
    Num_UEs = round(Total_Pop / people_per_ue);
    fprintf('Total Population in Area: %.0f\n', Total_Pop);
    fprintf('Generating %d UEs (1 UE per %d people)...\n', Num_UEs, people_per_ue);
    
    % 6. Create the Weighted Probability Distribution
    % Find the linear indices of all pixels that have actual people in them
    valid_pixels = find(local_data > 0);
    
    % Extract the population numbers for just those valid pixels to use as weights
    pixel_weights = local_data(valid_pixels); % weight based on population

    % 7. Sample the pixels! 
    % We use 'randsample' with replacement, weighted by the population.
    % If a pixel has 30,000 people, it will get picked roughly 100 times.
    sampled_pixel_indices = randsample(valid_pixels, Num_UEs, true, pixel_weights);
    
    % 8. Convert the 1D linear indices back to 2D matrix rows/cols
    [rows, cols] = ind2sub(size(local_data), sampled_pixel_indices);
    
    % 9. Convert the Rows/Cols back to real-world Latitudes and Longitudes
    [UE_lats_flat, UE_lons_flat] = intrinsicToGeographic(local_R, cols, rows);
    
    % 10. (Optional) Add spatial Jitter 
    % Because pixels are ~1km wide, multiple UEs in the same pixel will spawn 
    % on the exact same coordinate. We add a tiny random offset (±0.004 deg) 
    % so they scatter naturally within their 1km pixel bounds.
    jitter_deg = 0.004*3; 
    UE_lats_flat = UE_lats_flat + (rand(size(UE_lats_flat)) - 0.5) * jitter_deg;
    UE_lons_flat = UE_lons_flat + (rand(size(UE_lons_flat)) - 0.5) * jitter_deg;

    fprintf('Successfully generated geographic coordinates for %d UEs.\n', length(UE_lats_flat));
end