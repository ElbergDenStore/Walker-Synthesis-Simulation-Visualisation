function sir_lin_mat = Interference_calc_2D(el_mat, az_mat, BeamGrid, Cfg)
    NumUEs = size(az_mat, 1);
    nT = size(az_mat, 2);
    sir_lin_mat = nan(NumUEs, nT);

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
    
    served_mask = ~isnan(u_ues) & ~isnan(v_ues);
    
    for i = 1:NumUEs
        valid_t = find(served_mask(i,:));
        if isempty(valid_t), continue; end
        
        u = u_ues(i, valid_t)';
        v = v_ues(i, valid_t)';
        dist2_to_beams = (u - BeamGrid.b_u').^2 + (v - BeamGrid.b_v').^2;
        [min_dist2, mb_idx_local] = min(dist2_to_beams, [], 2);
        
        inside_main_beam = min_dist2 <= BeamGrid.r_beam^2;
        if ~any(inside_main_beam), continue; end
        
        valid_t = valid_t(inside_main_beam);
        u = u(inside_main_beam);
        v = v(inside_main_beam);
        min_dist2 = min_dist2(inside_main_beam);
        mb_idx_local = mb_idx_local(inside_main_beam);
        
        sig_lin = gain_from_dist2(min_dist2);
        
        nbs_matrix = BeamGrid.neighbor_idx(mb_idx_local, :); 
        valid_nbs_mask = ~isnan(nbs_matrix);
        safe_nbs_matrix = nbs_matrix;
        safe_nbs_matrix(~valid_nbs_mask) = 1;
        
        b_u_nbs = BeamGrid.b_u(safe_nbs_matrix); 
        b_v_nbs = BeamGrid.b_v(safe_nbs_matrix);
        
        d2_nb = (u - b_u_nbs).^2 + (v - b_v_nbs).^2;
        int_gains = gain_from_dist2(d2_nb);
        int_gains(~valid_nbs_mask) = 0;
        
        int_lin = sum(int_gains, 2) * Cfg.RU; 
        
        sir_lin_mat(i, valid_t) = sig_lin ./ max(int_lin, eps);
    end
end