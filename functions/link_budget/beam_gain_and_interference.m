function [carrier_density_dBmHz_mat, serving_beam_idx_mat, serving_beam_signal_density_dBmHz_mat, interference_density_dBmHz_mat] = beam_gain_and_interference(el_mat, az_mat, BeamGrid, Cfg)
% BEAM_GAIN_AND_INTERFERENCE  Serving-beam carrier and co-channel interference.
%   For each UE (rows) and time step (columns) in the elevation/azimuth matrices
%   EL_MAT/AZ_MAT, selects the serving beam from BEAMGRID and returns the carrier
%   and co-channel interference power spectral densities (dBm/Hz), the serving
%   beam index, and its signal density.  NaN inputs (no link) are returned as NaN.
    NumUEs = size(az_mat, 1);
    nT = size(az_mat, 2);
    
    % Initialize outputs
    carrier_density_dBmHz_mat = nan(NumUEs, nT);
    serving_beam_idx_mat      = nan(NumUEs, nT);
    serving_beam_signal_density_dBmHz_mat = nan(NumUEs, nT);
    interference_density_dBmHz_mat = nan(NumUEs, nT);

    Re = 6378.137e3;   
    
    % Step 1: Create a validity mask to only compute for active links (ignores NaNs)
    valid_mask = ~isnan(el_mat) & ~isnan(az_mat);
    valid_linear_idx = find(valid_mask); % Maps the flat 1D vectors back to the 2D matrix
    
    if isempty(valid_linear_idx)
        return; % Nothing to calculate
    end

    % Extract only valid points as flat column vectors
    el_flat = el_mat(valid_linear_idx);
    az_flat = az_mat(valid_linear_idx);
    
    el_flat = el_flat(:);
    az_flat = az_flat(:);
    
    % Step 2: Coordinate Transformation (Flat vectors)
    eta_flat = asind((Re / (Re + Cfg.Orbit_height)) .* cosd(el_flat));
    u_flat = sind(eta_flat) .* cosd(az_flat + 180);
    v_flat = sind(eta_flat) .* sind(az_flat + 180);

    % Ensure they are column vectors for knnsearch
    u_flat = u_flat(:);
    v_flat = v_flat(:);

    % Step 3: KD-Tree Nearest Neighbor Search (Serving Beam)
    % Ensure BeamGrid coordinates are column vectors
    beam_coords = [BeamGrid.u_center(:), BeamGrid.v_center(:)];
    query_coords = [u_flat, v_flat];

    [mb_idx_flat, ~] = knnsearch(beam_coords, query_coords);

    % Note: I removed the `inside_footprint` mask. The Array Factor naturally 
    % rolls off to massive losses outside the beam, making a hard mask unnecessary 
    % and saving computational overhead.

    % Step 4: Main Signal Gain
    du_serving = u_flat - BeamGrid.u_center(mb_idx_flat).';
    dv_serving = v_flat - BeamGrid.v_center(mb_idx_flat).';
    
    sig_loss_dB_flat = calculate_total_gain_dB(u_flat, v_flat, BeamGrid.u_center(mb_idx_flat).', BeamGrid.v_center(mb_idx_flat).', BeamGrid);
    
    % Center EIRP must be reshaped appropriately 
    center_eirp_flat = BeamGrid.BeamCenter_EIRP_dBmHz(mb_idx_flat).';
    sig_density_dBmHz_flat = center_eirp_flat + sig_loss_dB_flat; % loss is naturally negative

    % If ignoring interference, fill matrices and return early
    if isfield(Cfg, 'Ignore_Interference') && Cfg.Ignore_Interference
        carrier_density_dBmHz_mat(valid_linear_idx) = sig_density_dBmHz_flat;
        serving_beam_idx_mat(valid_linear_idx) = mb_idx_flat;
        serving_beam_signal_density_dBmHz_mat(valid_linear_idx) = sig_density_dBmHz_flat;
        interference_density_dBmHz_mat(valid_linear_idx) = -Inf;
        return;
    end

    % Step 5: Vectorized Neighbor Interference
    % Grab all neighbors for all valid data points simultaneously.
    % Assuming BeamGrid.neighbor_idx is [NumBeams x NumNeighbors]
    nbs_matrix = BeamGrid.neighbor_idx(mb_idx_flat, :); 
    
    valid_nbs_mask = ~isnan(nbs_matrix) & nbs_matrix > 0;
    safe_nbs_matrix = nbs_matrix;
    safe_nbs_matrix(~valid_nbs_mask) = 1; % Dummy index to prevent out-of-bounds crash

    % Pull neighbor coordinates and EIRP. Shape will be [N_valid x NumNeighbors]
    b_u_nbs = BeamGrid.u_center(safe_nbs_matrix); 
    b_v_nbs = BeamGrid.v_center(safe_nbs_matrix);
    nb_center_density = BeamGrid.BeamCenter_EIRP_dBmHz(safe_nbs_matrix);

    % Implicit expansion handles [N x 1] minus [N x NumNeighbors] cleanly
    du_nbs = u_flat - b_u_nbs;
    dv_nbs = v_flat - b_v_nbs;

    % Run the Array Factor pattern model on the entire N x NumNeighbors matrix
    int_loss_dB = calculate_total_gain_dB(u_flat, v_flat, b_u_nbs, b_v_nbs, BeamGrid);
    int_loss_dB(~valid_nbs_mask) = -Inf; % Erase fake neighbors by giving them infinite loss

    % Sum horizontally (dimension 2) and apply Resource Utilization (RU)
    int_contrib_mWHz = 10.^((nb_center_density + int_loss_dB) / 10);
    int_contrib_mWHz(~valid_nbs_mask) = 0; % Double check erasure in linear domain
    
    int_density_mWHz_flat = sum(int_contrib_mWHz, 2) * Cfg.RU;
    int_density_dBmHz_flat = 10 * log10(max(int_density_mWHz_flat, eps));

    % Step 6: Re-inflate back into the 2D Matrix using Linear Indexing
    carrier_density_dBmHz_mat(valid_linear_idx) = sig_density_dBmHz_flat;
    serving_beam_idx_mat(valid_linear_idx)      = mb_idx_flat;
    serving_beam_signal_density_dBmHz_mat(valid_linear_idx) = sig_density_dBmHz_flat;
    interference_density_dBmHz_mat(valid_linear_idx) = int_density_dBmHz_flat;
end

% --- Local Helper Function ---
function loss_dB = calculate_total_gain_dB(u_ue, v_ue, u_center, v_center, BeamGrid)
    % Array Factor (AF) - Same for both
    du = u_ue - u_center;
    dv = v_ue - v_center;
    
    du(du == 0) = eps; 
    dv(dv == 0) = eps;

    AF_u = sin(BeamGrid.Nu * (pi/2) .* du) ./ (BeamGrid.Nu .* sin((pi/2) .* du));
    AF_v = sin(BeamGrid.Nv * (pi/2) .* dv) ./ (BeamGrid.Nv .* sin((pi/2) .* dv));
    AF_loss_dB = 20 * log10(abs(AF_u .* AF_v));

    % Safely check for Element Factor defaults
    if ~isfield(BeamGrid, 'Cos_exponent')
        BeamGrid.Cos_exponent = 1.5; % Default flat patch
    end
    if ~isfield(BeamGrid, 'OneWeb')
        BeamGrid.OneWeb = false; % Default to Flat Phased Array
    end

    % Element Factor (EF) Routing
    if BeamGrid.OneWeb
        % === ONEWEB MODE ===
        % Element face is mechanically tilted to point directly at the beam center.
        % Roll-off is based strictly on the distance from the beam center (du, dv).
        rho_sq = du.^2 + dv.^2;
        rho_sq(rho_sq > 1) = 1; 
        
        cos_theta = sqrt(1 - rho_sq);
        EF_loss_dB = 10 * BeamGrid.Cos_exponent * log10(max(cos_theta, eps));
        
        loss_dB = AF_loss_dB + EF_loss_dB;
        
    else
        % === FLAT PHASED ARRAY MODE (Default) ===
        % Element face points at Nadir (u=0, v=0). We must normalize the loss 
        % so that the beam center is strictly 0 dB.
        
        % A. UE Element Loss (absolute from Nadir)
        rho_sq_ue = u_ue.^2 + v_ue.^2;
        rho_sq_ue(rho_sq_ue > 1) = 1;
        cos_theta_ue = sqrt(1 - rho_sq_ue);
        EF_ue_dB = 10 * BeamGrid.Cos_exponent * log10(max(cos_theta_ue, eps));
        
        % B. Center Element Loss (absolute from Nadir)
        rho_sq_center = u_center.^2 + v_center.^2;
        rho_sq_center(rho_sq_center > 1) = 1;
        cos_theta_center = sqrt(1 - rho_sq_center);
        EF_center_dB = 10 * BeamGrid.Cos_exponent * log10(max(cos_theta_center, eps));
        
        % C. Relative EF Loss
        relative_EF_loss_dB = EF_ue_dB - EF_center_dB;
        
        loss_dB = AF_loss_dB + relative_EF_loss_dB;
    end
end