function Link_2D = link_calc_matrix(el_mat, az_mat, range_mat, lat_vec, lon_vec, link_cfg, general_config)
    % Inputs are now NumUEs x nT matrices
    Link_2D.Frequency = link_cfg.f;
    Link_2D.Bandwidth = link_cfg.B;
    
    % FSPL
    c = physconst('LightSpeed');
    lambda = c / link_cfg.f;
    Link_2D.FSPL = 10 * log10(((4 * pi * range_mat) / lambda).^2);
    
    % Steering Loss
    Link_2D.Rx_steering_loss = Steering_loss_calc_2D(el_mat, link_cfg.Rx_type);
    
    % Absorption
    [Link_2D.Absorption_At, Link_2D.T_antenna] = Absorption_calc_2D(link_cfg.f, el_mat, lat_vec, lon_vec, general_config);
    
    % 4. Total Loss
    Link_2D.Total_loss = Link_2D.FSPL + Link_2D.Absorption_At + Link_2D.Rx_steering_loss;
    
    
    [carrier_density_dBmHz_mat, serving_beam_idx_mat, serving_beam_signal_density_dBmHz_mat, interference_density_dBmHz_mat] = beam_gain_and_interference(el_mat, az_mat, link_cfg.BeamGrid, general_config);
    
    area = 10*log10(4*pi.*range_mat.^2);
    Link_2D.PFD_W_MHz = carrier_density_dBmHz_mat - 30 - area - 10*log10(1e6);

    F_lin = 10^(link_cfg.NF/10);
    T_rx = (F_lin - 1) * 290;
    T_sys = Link_2D.T_antenna + T_rx;
    Link_2D.Noise_density_dBmHz = (10*log10(T_sys) + 10*log10(1.38e-23)) + 30;

    Link_2D.Carrier_density_dBmHz = carrier_density_dBmHz_mat + link_cfg.G_rx - Link_2D.Total_loss;
    Link_2D.Interference_density_dBmHz = interference_density_dBmHz_mat + link_cfg.G_rx - Link_2D.Total_loss;

    Link_2D.Rx_Power = Link_2D.Carrier_density_dBmHz;
    Link_2D.P_noise = Link_2D.Noise_density_dBmHz;
    Link_2D.interference_lin = 10.^(Link_2D.Interference_density_dBmHz / 10);
    Link_2D.serving_beam_idx = serving_beam_idx_mat;
    Link_2D.serving_beam_signal_lin = serving_beam_signal_density_dBmHz_mat;

    Link_2D.SNR = Link_2D.Carrier_density_dBmHz - Link_2D.Noise_density_dBmHz;
    Link_2D.SIR = Link_2D.Carrier_density_dBmHz - Link_2D.Interference_density_dBmHz;
    Link_2D.SINR = calculate_sinr_from_density(Link_2D.Carrier_density_dBmHz, Link_2D.Interference_density_dBmHz, Link_2D.Noise_density_dBmHz);

    if isfield(general_config, 'Ignore_Interference') && general_config.Ignore_Interference
        Link_2D.SIR = nan(size(el_mat));
        Link_2D.SINR = Link_2D.SNR;
        Link_2D.interference_lin = zeros(size(el_mat));
        Link_2D.Interference_density_dBmHz = -Inf(size(el_mat));
    end
end

% function BeamGrid = get_satellite_beam_grid(general_config)
%     if ~isfield(general_config, 'Satellite_antenna') || ~isfield(general_config.Satellite_antenna, 'BeamGrid') || isempty(general_config.Satellite_antenna.BeamGrid)
%         error('link_calc_matrix requires general_config.Satellite_antenna.BeamGrid to be populated before link calculation.');
%     end

%     BeamGrid = general_config.Satellite_antenna.BeamGrid;
% end

function SINR_dB = calculate_sinr_from_density(carrier_density_dBmHz, interference_density_dBmHz, noise_density_dBmHz)
    carrier_mWHz = 10.^(carrier_density_dBmHz / 10);
    interference_mWHz = 10.^(interference_density_dBmHz / 10);
    noise_mWHz = 10.^(noise_density_dBmHz / 10);

    SINR_lin = carrier_mWHz ./ (interference_mWHz + noise_mWHz);
    SINR_dB = 10 * log10(SINR_lin);
end

% function [adjusted_eirp_dBmHz, tx_steering_loss_dB] = calculate_beam_density_and_tx_steering(BeamGrid, serving_beam_idx)
%     if ~isfield(BeamGrid, 'BeamCenter_EIRP_dBmHz') || isempty(BeamGrid.BeamCenter_EIRP_dBmHz)
%         error('Satellite beam grid must contain BeamCenter_EIRP_dBmHz. Populate general_config.Satellite_antenna.BeamGrid during scenario setup.');
%     end

%     adjusted_eirp_dBmHz = nan(size(serving_beam_idx));
%     tx_steering_loss_dB = nan(size(serving_beam_idx));

%     valid_mask = isfinite(serving_beam_idx) & serving_beam_idx > 0;
%     if ~any(valid_mask)
%         return;
%     end

%     adjusted_eirp_dBmHz(valid_mask) = BeamGrid.BeamCenter_EIRP_dBmHz(serving_beam_idx(valid_mask));
%     tx_steering_loss_dB(valid_mask) = BeamGrid.TxSteeringLoss_dB(serving_beam_idx(valid_mask));
% end

% --- Subfunctions adjusted for matrices ---
function loss_rx = Steering_loss_calc_2D(el_mat, rx_type)
    theta_rx = 90 - el_mat;

    loss_rx = zeros(size(el_mat));
    if contains(rx_type, 'array')
        loss_rx = -10 * log10(cosd(theta_rx).^1.5);
    end
end

function [At_mat, Tsky_mat] = Absorption_calc_2D(f, el_mat, lat_vec, lon_vec, Cfg)
    % Fast Mode / Fallback: Instantaneous matrix fill
    static_loss = 0.5;
    if isfield(Cfg, 'Simple_Atmospheric_Loss_dB'), static_loss = Cfg.Simple_Atmospheric_Loss_dB; end

    if isfield(Cfg, 'Use_P618') && Cfg.Use_P618 == false
        At_mat = zeros(size(el_mat)) + static_loss;
        Tsky_mat = zeros(size(el_mat)) + 290;
        return;
    end

    % Preferred precision path: direct LUT interpolation only!
    freq_GHz = f / 1e9;
    lut_file_target = sprintf('p618_%.1f.mat', freq_GHz);
    
    if isfield(Cfg, 'P618_LUT_File') && strlength(string(Cfg.P618_LUT_File)) > 0
        lut_file_target = char(string(Cfg.P618_LUT_File));
    end

    lut_path = which(lut_file_target);

    % If exact match not found, scan functions/data for best frequency match within 20%
    if isempty(lut_path)
        script_dir = fileparts(mfilename('fullpath'));
        data_dir = fullfile(script_dir, 'data');
        candidate_luts = dir(fullfile(data_dir, 'p618_*.mat'));
        best_err = Inf;
        best_lut_freq_hz = NaN;
        for ci = 1:numel(candidate_luts)
            cpath = fullfile(data_dir, candidate_luts(ci).name);
            % Check if file contains a 'LUT' variable before loading (avoids spurious warnings)
            file_info = whos('-file', cpath, 'LUT');
            if isempty(file_info)
                continue;
            end
            try
                Sc = load(cpath, 'LUT');
            catch
                continue;
            end
            if ~isfield(Sc, 'LUT') || ~isfield(Sc.LUT, 'frequency_hz')
                continue;
            end
            err = abs(Sc.LUT.frequency_hz - f) / f;
            if err < best_err
                best_err = err;
                lut_path = cpath;
                best_lut_freq_hz = Sc.LUT.frequency_hz;
            end
        end
        if ~isempty(lut_path) && best_err > 0.20
            warning('LinkCalc:P618FrequencyMismatch', 'Closest LUT is %.1f GHz but requested %.1f GHz (%.0f%% off). Falling back to simple attenuation.', ...
                best_lut_freq_hz/1e9, freq_GHz, best_err*100);
            lut_path = '';
        elseif ~isempty(lut_path)
            fprintf('P618: Using %.1f GHz LUT for %.1f GHz request (%.1f%% freq difference).\n', ...
                best_lut_freq_hz/1e9, freq_GHz, best_err*100);
        end
    end

    lut_valid = false;
    if ~isempty(lut_path) && isfile(lut_path)
        S = load(lut_path, 'LUT');
        if isfield(S, 'LUT')
            [At_mat, Tsky_mat, ok] = interpolate_p618_lut_simple(S.LUT, el_mat, lat_vec, lon_vec);
            if ok
                lut_valid = true;
            end
        end
    else
        warning('LinkCalc:P618LUTNotFound', 'No suitable LUT found for %.1f GHz. Falling back to simple attenuation.', freq_GHz);
    end

    if ~lut_valid
        At_mat = zeros(size(el_mat)) + static_loss;
        Tsky_mat = zeros(size(el_mat)) + 290;
    end
end

function [At_mat, Tsky_mat, ok] = interpolate_p618_lut_simple(LUT, el_mat, lat_vec, lon_vec)
    ok = false;
    At_mat = [];
    Tsky_mat = [];

    required = {'lat_deg', 'lon_deg', 'el_deg', 'At_dB', 'Tsky_K'};
    for i = 1:numel(required)
        if ~isfield(LUT, required{i})
            return;
        end
    end

    lat_min = min(LUT.lat_deg);
    lat_max = max(LUT.lat_deg);
    lon_min = min(LUT.lon_deg);
    lon_max = max(LUT.lon_deg);

    req_lat_min = min(lat_vec(:));
    req_lat_max = max(lat_vec(:));
    req_lon_min = min(lon_vec(:));
    req_lon_max = max(lon_vec(:));

    % Produce warning and fallback if out of bounds
    if req_lat_min < lat_min - 0.1 || req_lat_max > lat_max + 0.1 || ...
       req_lon_min < lon_min - 0.1 || req_lon_max > lon_max + 0.1
        warning('LinkCalc:P618LUTOutOfBounds', ...
                'Requested UE coordinates (Lat: %.1f to %.1f, Lon: %.1f to %.1f) are outside the LUT generated range (Lat: %.1f to %.1f, Lon: %.1f to %.1f). Falling back to simple attenuation.', ...
                req_lat_min, req_lat_max, req_lon_min, req_lon_max, lat_min, lat_max, lon_min, lon_max);
        return;
    end

    n_ues = numel(lat_vec);
    n_t = size(el_mat, 2);

    el_min = min(LUT.el_deg);
    el_max = max(LUT.el_deg);

    % Safely clamp coordinates slightly to avoid griddedInterpolant NaN at the exact edges
    lat_q = repmat(min(max(lat_vec(:), lat_min), lat_max), 1, n_t);
    lon_q = repmat(min(max(lon_vec(:), lon_min), lon_max), 1, n_t);
    el_q = min(max(el_mat, el_min), el_max);

    try
        F_at = griddedInterpolant({double(LUT.lat_deg), double(LUT.lon_deg), double(LUT.el_deg)}, ...
            double(LUT.At_dB), 'linear', 'nearest');
        F_tsky = griddedInterpolant({double(LUT.lat_deg), double(LUT.lon_deg), double(LUT.el_deg)}, ...
            double(LUT.Tsky_K), 'linear', 'nearest');

        At_mat = F_at(double(lat_q), double(lon_q), double(el_q));
        Tsky_mat = F_tsky(double(lat_q), double(lon_q), double(el_q));

        At_mat = reshape(At_mat, [n_ues, n_t]);
        Tsky_mat = reshape(Tsky_mat, [n_ues, n_t]);
        ok = true;
    catch
        ok = false;
    end
end