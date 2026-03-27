function [sir_lin_mat, mb_idx_mat] = Interference_calc_2D(el_mat, az_mat, BeamGrid, Cfg)
    NumUEs = size(az_mat, 1);
    nT = size(az_mat, 2);
    sir_lin_mat = nan(NumUEs, nT);
    mb_idx_mat  = nan(NumUEs, nT);

    Re = 6378.14e3;   
    
    eta_mat = asind((Re / (Re + Cfg.Orbit_height)) * cosd(el_mat));
    
    % Convert to u and v coordinates for the entire matrix
    u_ues = sind(eta_mat) .* cosd(az_mat + 180);
    v_ues = sind(eta_mat) .* sind(az_mat + 180);
    
    N_elements = round(101.2 / BeamGrid.Beamwidth_deg); 
    d_spacing = BeamGrid.lambda / 2;
    phase_term = @(d2) (pi * d_spacing / BeamGrid.lambda) * max(sqrt(max(0, d2)), eps);
    array_factor = @(d2) sin(N_elements * phase_term(d2)) ./ (N_elements * sin(phase_term(d2)));
    gain_from_dist2 = @(d2) array_factor(d2).^2;

    %%%%%%%% To avoid the for loop, a knnsearch approach is made. A
    %%%%%%%% "vector" version would create a beams x ues x time matrix that
    %%%%%%%% would be crazy...

    %%%%%%%% for 108 UEs and 12h, a 4x improvement of computational speed. the bigger, the larger improvement. 
    
    valid_mask = ~isnan(u_ues) & ~isnan(v_ues);
    valid_linear_idx = find(valid_mask); % Store exact original matrix locations

    u_flat = u_ues(valid_linear_idx);
    v_flat = v_ues(valid_linear_idx);

    % Step B: KD-Tree Nearest Neighbor Search (The Magic Trick)
    % This finds the closest beam for millions of points without blowing up RAM
    beam_coords = [BeamGrid.b_u, BeamGrid.b_v];
    query_coords = [u_flat, v_flat];

    [mb_idx_flat, min_dist_flat] = knnsearch(beam_coords, query_coords);
    min_dist2_flat = min_dist_flat.^2;

    % Step C: Filter out UEs that fall outside the main beam footprint
    inside_mask = min_dist2_flat <= BeamGrid.r_beam^2;

    % Keep only the data where the UE is actually served by a beam
    final_linear_idx = valid_linear_idx(inside_mask); % Final map back to the 2D matrix
    u_final          = u_flat(inside_mask);
    v_final          = v_flat(inside_mask);
    mb_idx_final     = mb_idx_flat(inside_mask);
    min_dist2_final  = min_dist2_flat(inside_mask);

    % Step D: Main Signal Gain
    sig_lin_final = gain_from_dist2(min_dist2_final);

    % Step E: Vectorized Neighbor Interference
    % Grab all 6 neighbors for all millions of valid data points simultaneously
    nbs_matrix = BeamGrid.neighbor_idx(mb_idx_final, :); 
    valid_nbs_mask = ~isnan(nbs_matrix);
    safe_nbs_matrix = nbs_matrix;
    safe_nbs_matrix(~valid_nbs_mask) = 1;

    b_u_nbs = BeamGrid.b_u(safe_nbs_matrix); 
    b_v_nbs = BeamGrid.b_v(safe_nbs_matrix);

    % Distance squared to all 6 neighbors (MATLAB automatically expands this N x 6)
    d2_nb = (u_final - b_u_nbs).^2 + (v_final - b_v_nbs).^2;

    % Run the phased array math on the entire N x 6 neighbor matrix
    int_gains = gain_from_dist2(d2_nb);
    int_gains(~valid_nbs_mask) = 0; % Erase fake neighbors

    % Sum horizontally and apply Resource Utilization
    int_lin_final = sum(int_gains, 2) * Cfg.RU; 

    % Step F: Final SIR Calculation
    sir_final = sig_lin_final ./ max(int_lin_final, eps);

    % Step G: Re-inflate back into the 2D Matrix using Linear Indexing!
    sir_lin_mat(final_linear_idx) = sir_final;

    mb_idx_mat(final_linear_idx)  = mb_idx_final;
end