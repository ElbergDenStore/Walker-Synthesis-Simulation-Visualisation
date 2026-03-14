function [UE_lats_flat, UE_lons_flat] = generate_equal_ish_area_UEs(Lat_vec, Lon_vec)
    lat_min = min(Lat_vec);
    lat_max = max(Lat_vec);
    lon_min = min(Lon_vec);
    lon_max = max(Lon_vec);
        
     
        
    % 1. Calculate how many latitude rows we need
    lat_dist_km = (lat_max - lat_min) * 111.32; % ~111.32 km per deg of latitude
    num_lats = length(Lat_vec);
    % Define your target physical spacing between UEs
    target_spacing_km = lat_dist_km/num_lats;
    
    % num_lats = max(2, round(lat_dist_km / target_spacing_km) + 1); 
    Lat_vec = linspace(lat_min, lat_max, num_lats); % Guarantees lat limits are hit
    
    UE_lats_flat = [];
    UE_lons_flat = [];
    
    % 2. Build the dynamic longitude rows
    for i = 1:length(Lat_vec)
        current_lat = Lat_vec(i);
        
        % Calculate physical width of this specific longitude span
        lon_span_km = (lon_max - lon_min) * 111.32 * cosd(current_lat);
        
        % Calculate how many points we need to maintain the target spacing
        num_lons = max(2, round(lon_span_km / target_spacing_km) + 1);
        
        % Guarantees -60 and 30 are hit on every row
        current_lons = linspace(lon_min, lon_max, num_lons); 
        
        % Append to our flat lists
        UE_lats_flat = [UE_lats_flat, repmat(current_lat, 1, num_lons)];
        UE_lons_flat = [UE_lons_flat, current_lons];
    end
end