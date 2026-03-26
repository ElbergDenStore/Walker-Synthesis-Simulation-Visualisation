function BeamGrid = calculate_Beams(f, G_tx, orbit_height, Min_Elev_deg, FRF)
    % Default parameters if not explicitly provided
    if nargin < 4, Min_Elev_deg = 20; end
    if nargin < 5, FRF = 1; end
    
    % 1. Physics & Constants
    Re = 6371; % Earth radius in km
    h = orbit_height / 1000; % Convert meters to km
    c = physconst('LightSpeed'); 
    lambda = c / f;
    
    % 2. Beam geometry
    Beamwidth_deg = sqrt(32400 ./ (10.^(G_tx/10)));
    r_beam = sind(Beamwidth_deg / 2);
    eta_max = asind((Re / (Re + h)) * cosd(Min_Elev_deg));
    du = sind(Beamwidth_deg) / 2;
    rings = ceil(sind(eta_max) / du) + 2;
    
    % 3. Hex Grid Generation
    [Q, R] = meshgrid(-rings:rings, -rings:rings);
    Q = Q(:); 
    R = R(:);
    
    % Calculate u and v for the entire grid simultaneously
    U = du * sqrt(3) * (Q + R/2); 
    V = du * 1.5 * R;
    
    % Create a logical mask for beams inside the field of view
    valid_mask = (U.^2 + V.^2) <= sind(eta_max)^2;
    
    % Filter down to only the valid beams
    b_u = U(valid_mask);
    b_v = V(valid_mask);
    b_q = Q(valid_mask);
    b_r = R(valid_mask);
    
    num_beams = length(b_u);
    
    % 4. Neighbor Mapping based on Frequency Reuse Factor (FRF)
    switch FRF
        case 1
            neighbor_offsets = [1 0; -1 0; 0 1; 0 -1; 1 -1; -1 1];
        case 3
            neighbor_offsets = [2 -1; 1 1; -1 2; -2 1; -1 -1; 1 -2];
        case 4
            neighbor_offsets = [2 0; -2 0; 0 2; 0 -2; 2 -2; -2 2];
        otherwise
            error('Unsupported FRF. Please use 1, 3, or 4.');
    end
    
    valid_coords = [b_q, b_r];
    neighbor_idx = nan(num_beams, 6);
    
    % Instead of looping over every beam, loop over the 6 offsets
    for k = 1:6
        % Shift EVERY beam by the current offset simultaneously
        target_coords = valid_coords + neighbor_offsets(k, :);
        
        % ismember instantly finds which of our target coords actually exist in our grid
        % 'found' is a logical array, 'loc' is the exact row index of the neighbor
        [found, loc] = ismember(target_coords, valid_coords, 'rows');
        
        % Assign the found locations directly into our neighbor matrix
        temp_idx = nan(num_beams, 1);
        temp_idx(found) = loc(found);
        neighbor_idx(:, k) = temp_idx;
    end
    
    % 5. Bundle everything into the output struct
    BeamGrid.b_u = b_u;
    BeamGrid.b_v = b_v;
    BeamGrid.b_q = b_q;
    BeamGrid.b_r = b_r;
    BeamGrid.neighbor_idx = neighbor_idx;
    BeamGrid.r_beam = r_beam;
    BeamGrid.lambda = lambda;
    BeamGrid.Beamwidth_deg = Beamwidth_deg;
    BeamGrid.num_beams = num_beams;
    BeamGrid.FRF = FRF;
end