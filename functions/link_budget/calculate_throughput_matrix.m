function Throughput = calculate_throughput_matrix(SINR_dB, serving_beam_idx_mat, serving_sat_id_mat, BeamGrid, bandwidth_Hz, share_bandwidth, modified_shannon)
    %% TODO, change the name so matrix is not the end of all the functions. the variables is cool to have, but not the functions
    Throughput = zeros(size(SINR_dB));
    valid_idx = ~isnan(SINR_dB);

    SINR_lin_valid = 10.^(SINR_dB(valid_idx) / 10);

    if modified_shannon
        %% P. Mogensen et al. "LTE Capacity Compared to the Shannon Bound", 2007 IEEE 65th Vehicular Technology Conference - VTC2007-Spring, Dublin, Ireland, 2007, pp. 1234-1238, doi: 10.1109/VETECS.2007.260.
        BW_eff = 0.56;
        eta = 1;
        SNR_eff_dB = 2;
        SNR_eff = 10^(SNR_eff_dB / 10);
    else 
        BW_eff  = 1; %% Normal shannon capacity
        eta     = 1;
        SNR_eff = 1;
    end
    raw_throughput = bandwidth_Hz * BW_eff * eta .* log2(1 + SINR_lin_valid ./ SNR_eff);
    
    if share_bandwidth
        [num_ues, n_t] = size(serving_beam_idx_mat);
        time_matrix = repmat(1:n_t, num_ues, 1);

        valid_mask = ~isnan(serving_beam_idx_mat);
        beam_valid = serving_beam_idx_mat(valid_mask);
        time_valid = time_matrix(valid_mask);
        sat_valid = serving_sat_id_mat(valid_mask);

        if isfield(BeamGrid, 'Total_beams') % Using new pre-computed standard grids if available
            num_beams_for_id = BeamGrid.Total_beams;
        elseif isfield(BeamGrid, 'num_beams')
            num_beams_for_id = BeamGrid.num_beams;
        else
            num_beams_for_id = 10000; % Safe arbitrary large spacing
        end

        % We map (beam_idx, time_idx, sat_idx) to ensure we don't accidentally class UEs under different satellites as sharing the same beam.
        % Adding sat_idx multiplied by a very large spacer ensures zero collision across satellites.
        max_time_steps = n_t + 1;
        unique_beam_time_ids = beam_valid + (time_valid * num_beams_for_id) + (sat_valid * num_beams_for_id * max_time_steps);
        [~, ~, ic] = unique(unique_beam_time_ids);
        counts = accumarray(ic, 1);
        users_sharing_beam = counts(ic);

        users_per_beam_mat = ones(num_ues, n_t);
        users_per_beam_mat(valid_mask) = users_sharing_beam;
        Throughput(valid_idx) = raw_throughput ./ users_per_beam_mat(valid_idx);
    else
        Throughput(valid_idx) = raw_throughput;
    end
end