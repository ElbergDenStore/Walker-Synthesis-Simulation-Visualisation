function [UE_lats_flat, UE_lons_flat] = generate_random_UEs(lat_vec, lon_vec, num_UEs)
    latlim = [min(lat_vec), max(lat_vec)];
    lonlim = [min(lon_vec), max(lon_vec)];
    nLat = 1e3; nLon = 1e3; 
    [LonG, LatG] = meshgrid(linspace(lonlim(1), lonlim(2), nLon), linspace(latlim(1), latlim(2), nLat));

    sampled_pixel_indices = randsample(LatG, num_UEs);
    
    UE_lats_flat=LatG(sampled_pixel_indices);
    UE_lons_flat=LonG(sampled_pixel_indices);
end