function [UE_lats_flat, UE_lons_flat] = generate_equal_area_ues(lat_limits, lon_limits, num_ues)
% GENERATE_SPHERICAL_UES Creates uniformly spaced UEs accounting for true spherical
% geometry at high latitudes, while strictly guaranteeing the exact num_ues.
    lat_min = min(lat_limits);
    lat_max = max(lat_limits);
    lon_min = min(lon_limits);
    lon_max = max(lon_limits);
    num_ues = max(1, round(num_ues));

    %% 1. SPHERICAL GEOMETRY (The Shape)
    R = 6371; % Earth radius in km
    lat_rad = deg2rad([lat_min, lat_max]);
    lon_rad = deg2rad([lon_min, lon_max]);
    
    % Exact spherical area of the bounding box (km^2)
    A = R^2 * (lon_rad(2) - lon_rad(1)) * (sin(lat_rad(2)) - sin(lat_rad(1)));
    
    % Ideal physical spacing between UEs (km)
    d = sqrt(A / num_ues);
    
    % True North-South arc length (km)
    H = R * (lat_rad(2) - lat_rad(1));
    
    % Determine the optimal number of rows to maintain the spacing 'd'
    num_lat_rows = max(1, round(H / d));
    num_lat_rows = min(num_lat_rows, num_ues); % Safety check for low UE counts
    
    % Create the latitude rows (hitting the exact limits as preferred)
    lat_rows = linspace(lat_min, lat_max, num_lat_rows);

    %% 2. STRICT ACCOUNTING (The Exact Count)
    % Weight rows by the physical width of the Earth at that latitude
    row_weights = cosd(lat_rows);
    row_weights = max(row_weights, 0);
    if sum(row_weights) == 0
        row_weights = ones(size(row_weights));
    end
    
    % Calculate exact fractional points needed per row
    ideal_counts = (row_weights / sum(row_weights)) * num_ues;
    
    % Base allocation (round down to guarantee we don't overshoot)
    row_counts = floor(ideal_counts);
    
    % How many UEs are left over after rounding down?
    count_delta = num_ues - sum(row_counts);
    
    % Distribute the exact remainder to hit num_ues perfectly
    if count_delta > 0
        % Sort by which rows lost the most fractional value when rounding down
        [~, add_order] = sort(ideal_counts - row_counts, 'descend');
        for k = 1:count_delta
            % Distribute one UE at a time to the most deserving rows
            idx = add_order(mod(k - 1, numel(add_order)) + 1);
            row_counts(idx) = row_counts(idx) + 1;
        end
    end

    %% 3. BUILD THE FLAT ARRAYS
    UE_lats_flat = zeros(num_ues, 1);
    UE_lons_flat = zeros(num_ues, 1);
    write_idx = 1;

    for r = 1:numel(lat_rows)
        n_row = row_counts(r);
        
        % Skip if the accounting phase gave this row 0 points
        if n_row <= 0
            continue;
        end

        % Create longitudes for this row
        if n_row == 1
            lon_row = mean([lon_min, lon_max]);
        else
            lon_row = linspace(lon_min, lon_max, n_row);
        end

        % Write to the flat pre-allocated arrays
        idx_range = write_idx:(write_idx + n_row - 1);
        UE_lats_flat(idx_range) = lat_rows(r);
        UE_lons_flat(idx_range) = lon_row;
        
        write_idx = write_idx + n_row;
    end
end