function BeamGrid = calculate_hexagonal_beams(G_tx_dBi, f_Hz, orbit_height_m, Min_Elev_deg, target_EIRP_or_PFD, FRF)
%CALCULATE_HEXAGONAL_BEAMS Generate mathematical model of a flat phased array
% creating a hexagonal grid of spot beams across the earth footprint.
    if nargin < 4 || isempty(Min_Elev_deg), Min_Elev_deg = 20; end
    if nargin < 5 || isempty(target_EIRP_or_PFD), target_EIRP_or_PFD = -128; end % Treat as PFD target by default
    if nargin < 6 || isempty(FRF), FRF = 1; end

    % Array Physics
    BeamGrid.OneWeb = false; 
    BeamGrid.Cos_exponent = 1.5; 
    BeamGrid.Element_gain = 6; % dBi

    % Reverse-calculate the number of elements needed to hit the target Peak Gain
    % G_total = 10*log10(N) + G_element  => N = 10^((G - G_el)/10)
    total_elements = 10^((G_tx_dBi - BeamGrid.Element_gain) / 10);
    BeamGrid.Nu = max(round(sqrt(total_elements)), 1);
    BeamGrid.Nv = max(round(sqrt(total_elements)), 1);
    
    BeamGrid.Max_gain = 10*log10(BeamGrid.Nu * BeamGrid.Nv) + BeamGrid.Element_gain;

    % Theoretical 3dB beamwidth in u/v direction cosine space
    uv_3dB = 2 * (1.391 / (BeamGrid.Nu * (pi/2))); 
    r_beam = uv_3dB / 2; % radius

    % Hexagonal center layout mechanics
    Re_km = 6371; 
    h_km = orbit_height_m / 1000; 
    eta_max = asind((Re_km / (Re_km + h_km)) * cosd(Min_Elev_deg));

    % Generate hex rings up to edge of earth
    du = sqrt(3) * r_beam; % spacing for tight hexagonal packing
    rings = ceil(sind(eta_max) / du) + 1;

    [Q, R] = meshgrid(-rings:rings, -rings:rings);
    Q = Q(:); R = R(:);

    % Convert axial coords to cartesian U/V
    U = sqrt(3)*r_beam * (Q + R/2);
    V = 1.5*r_beam * R;

    % Clip to visible bounds
    valid_mask = (U.^2 + V.^2) <= sind(eta_max)^2;
    b_u = U(valid_mask);
    b_v = V(valid_mask);
    b_q = Q(valid_mask);
    b_r = R(valid_mask);
    beam_count = length(b_u);

    BeamGrid.num_beams = beam_count;
    BeamGrid.u_center = b_u.';
    BeamGrid.v_center = b_v.';

    % Frequency Reuse factor channel map
    if FRF == 1
        channel_id = ones(1, beam_count);
    elseif FRF == 3
        channel_id = mod(b_q + 2 * b_r, 3)' + 1;
    elseif FRF == 4
        channel_id = mod(b_q, 2)' + 2 * mod(b_r, 2)' + 1;
    else
        channel_id = mod(1:beam_count, FRF) + 1;
    end

    % Adjust power scaling 
    theta_rad = asin(sqrt(min(b_u.^2 + b_v.^2, 1)));
    R_sat = (Re_km + h_km) * 1000;
    slant_range_m = R_sat .* cos(theta_rad) - sqrt((Re_km*1000)^2 - (R_sat .* sin(theta_rad)).^2);
    area_spreading_dB = 10 * log10(4 * pi * (slant_range_m.^2));

    if target_EIRP_or_PFD < -50
        % Input was likely a PFD goal (-115 dBW/m2 etc), back calculate EIRP density required at bore
        BeamGrid.BeamCenter_EIRP_dBmHz = (target_EIRP_or_PFD + 30) + area_spreading_dB;
    else
        % Input was EIRP, just apply it
        BeamGrid.BeamCenter_EIRP_dBmHz = repmat(target_EIRP_or_PFD, 1, beam_count);
    end

    % Calculate valid interference neighbors matching the channel ID
    max_neighbors = 15; % Preallocate
    BeamGrid.neighbor_idx = nan(beam_count, max_neighbors);
    
    for i = 1:beam_count
        my_chan = channel_id(i);
        % Find everyone else on the same channel
        fellows = find(channel_id == my_chan);
        fellows(fellows == i) = []; % remove self
        
        % Sort them by physical u/v distance (closest interferers first)
        dists = sqrt((b_u(fellows) - b_u(i)).^2 + (b_v(fellows) - b_v(i)).^2);
        [~, sort_idx] = sort(dists, 'ascend');
        fellows = fellows(sort_idx);
        
        % Keep top neighbors
        n_keep = min(length(fellows), max_neighbors);
        BeamGrid.neighbor_idx(i, 1:n_keep) = fellows(1:n_keep);
    end
end