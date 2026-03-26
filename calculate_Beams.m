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
    % Preallocating to the maximum possible size and trimming is faster, 
    % but for a few hundred elements, dynamic growth is fine and clean.
    b_u = []; b_v = [];
    b_q = []; b_r = [];
    
    for q = -rings:rings
        for r = -rings:rings
            u_val = du * sqrt(3) * (q + r/2); 
            v_val = du * 1.5 * r;
            
            % Only keep beams that fall within the satellite's maximum field of view
            if (u_val^2 + v_val^2) <= sind(eta_max)^2
                b_u = [b_u; u_val]; 
                b_v = [b_v; v_val]; 
                b_q = [b_q; q];     
                b_r = [b_r; r];     
            end
        end
    end
    
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
    
    % Create a hashmap to instantly look up beam indices by their q,r coordinates
    beam_key_to_idx = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for b = 1:num_beams
        beam_key_to_idx(sprintf('%d_%d', b_q(b), b_r(b))) = b;
    end
    
    % Find the indices of all 6 neighbors for every beam
    neighbor_idx = nan(num_beams, 6);
    for b = 1:num_beams
        q0 = b_q(b);
        r0 = b_r(b);
        for k = 1:6
            qn = q0 + neighbor_offsets(k,1);
            rn = r0 + neighbor_offsets(k,2);
            key = sprintf('%d_%d', qn, rn);
            if isKey(beam_key_to_idx, key)
                neighbor_idx(b, k) = beam_key_to_idx(key);
            end
        end
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