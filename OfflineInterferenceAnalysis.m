
clearvars; close all; clc;

% %% 1. Configuration & Region Selection
% % Toggle between 'Nordjylland', 'Denmark', or 'Full'
% REGION = 'Full3000'; 

% switch REGION
%     case 'Nordjylland'
%         latlim = [56.5, 58.0];
%         lonlim = [8.0, 11.0];
%         people_per_ue = 300;
%         dataFile = 'Nordjylland300.mat';
%     case 'Denmark'
%         latlim = [54.5, 58.0];
%         lonlim = [8.0, 15.5];
%         people_per_ue = 3000;
%         dataFile = 'Denmark3000.mat';
%     case 'Full3000'
%         latlim = [54.5, 83.9]; 
%         lonlim = [-60, 30.0];
%         people_per_ue = 3000;
%         dataFile = 'Full3000.mat';
%     case 'Full30000'
%         latlim = [54.5, 83.9]; 
%         lonlim = [-60, 30.0];
%         people_per_ue = 30000;
%         dataFile = 'Full30000.mat';
% end

% FORCE_RERUN = false; 

% if exist(dataFile, 'file') && ~FORCE_RERUN 
%     fprintf('Loading %s data from %s...\n', REGION, dataFile);
%     load(dataFile);
% else
%     fprintf('Running new simulation for %s...\n', REGION);
%     constellation_type = "walkerDelta";
%     ue_grid_size = "medium";
%     duration = "short";
%     frequency = "fr1";
%     height_km = 1000;
%     Cfg = get_cfg(height_km);

    
%     Cfg.Accept_Flat_UE_array = true; Cfg.Equal_UE_area = false; 
%     [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, Total_Pop] = generate_population_based_UEs(latlim, lonlim, people_per_ue);

%     calc_link = true; 
%     Cfg.Use_P618 = false; %simple atmospheric loss
%     Cfg.Simple_Atmospheric_Loss_dB = 1;
%     plot_results = false; 
%     use_parallel = true;
%     metrics = coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link);
%     save(dataFile);
% end

% %%%%% VISUALISE UES
% % --- VISUAL VERIFICATION PLOT ---
% fprintf('Plotting UE distribution map...\n');

% % Create a clean, white-background figure window
% figure('Name',  'UE Distribution', 'Color', 'w', 'Position', [100, 100, 800, 800]);

% % Create geographic axes
% gx = geoaxes;

% % Scatter the UEs: Size 30, Red, Filled, with a slight black edge for visibility
% geoscatter(gx, Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, 30, 'r', 'filled', ...
%     'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% % Lock the map view specifically to Greenland's coordinates
% geolimits(gx, latlim, lonlim);

% % Add a basemap to show the terrain and ice sheet 
% % (Other good options: 'topographic', 'streets-light', or 'bluegreen')
% geobasemap(gx, 'satellite'); 

% % Add a dynamic title
% Cfg.NumUEs = length(Cfg.Flat_UE_array.Lats);
% title(sprintf(' UE Distribution\nTotal UEs: %d (1 per %d people)', ...
%     Cfg.NumUEs, round(Total_Pop/Cfg.NumUEs)), 'FontWeight', 'bold', 'FontSize', 14);

% pause(3*(Cfg.NumUEs/3000)) %to load UE distribution, 3000UEs take 3 seconds
% exportgraphics(gcf, 'screenshots/UEDistribution.png', 'Resolution', 300);


constellation_type = "walkerDelta";
ue_grid_size = "small";
duration = "short";
frequency = "fr1";
height_km = 1000;
Cfg = get_cfg(height_km);


calc_link = true; 
Cfg.Use_P618 = false; %simple atmospheric loss
Cfg.Simple_Atmospheric_Loss_dB = 1;
plot_results = false; 
use_parallel = false;
metrics = coverage_simulator_function(Cfg, plot_results, use_parallel, calc_link);
Cfg.NumUEs = size(metrics.Num_visible,1);




%% 3. Geometry & Beam Grid %%%%%%%%%%%% THESE ARE ALWAYS THE EXACT SAME %%%%%%%%%%%% IMAGINE THE SATELLITE ALWAYS ORIENTED THE SAME WAY
Re = 6371; h = Cfg.Orbit_height/1000; Min_Elev_deg = 20;
f = 20e9; c = 3e8; lambda = c/f; G = 40;
Beamwidth_deg = sqrt(32400./(10.^(G/10)));
r_beam = sind(Beamwidth_deg / 2);
eta_max = asind((Re / (Re + h)) * cosd(Min_Elev_deg));
du = sind(Beamwidth_deg) / 2;
rings = ceil(sind(eta_max) / du) + 2;
% Add arrays to store the hex grid logical coordinates
b_u = []; b_v = []; % preallocate to rings^2 or something and afterwards reduce
b_q = []; b_r = [];
for q = -rings:rings
    for r = -rings:rings
        u_val = du * sqrt(3) * (q + r/2); 
        v_val = du * 1.5 * r;
        if (u_val^2 + v_val^2) <= sind(eta_max)^2
            b_u = [b_u; u_val]; 
            b_v = [b_v; v_val];
            b_q = [b_q; q]; % Save q-index
            b_r = [b_r; r]; % Save r-index
        end
    end
end



num_beams = length(b_u)
%% 4. Vectorize UE Projections
nT = numel(metrics.SimData(1).Elevation_deg);
u_ues = nan(Cfg.NumUEs, nT);
v_ues = nan(Cfg.NumUEs, nT);
served_mask = false(Cfg.NumUEs, nT);

for i = 1:Cfg.NumUEs % make parfor soon or vector
    el = metrics.SimData(i).Elevation_deg;
    az = metrics.SimData(i).Azimuth_deg;
    served_mask(i,:) = ~isnan(el) & ~isnan(az);
    eta_vec = asind((Re / (Re + h)) * cosd(el));
    
    % nadir steering angles
    u_ues(i,:) = sind(eta_vec) .* cosd(az + 180);
    v_ues(i,:) = sind(eta_vec) .* sind(az + 180);
end

%% 5. Interference Model (Main Beam + 6 Hex Neighbors)
neighbor_offsets = [1 0; -1 0; 0 1; 0 -1; 1 -1; -1 1];
beam_key_to_idx = containers.Map('KeyType', 'char', 'ValueType', 'double');
for b = 1:num_beams
    beam_key_to_idx(sprintf('%d_%d', b_q(b), b_r(b))) = b;
end

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

gain_exp = 1.5;
% gain_from_dist2 = @(d2) max(cosd(asind(min(1, sqrt(max(0, d2))))), 0).^gain_exp; Complicated but safe version of the next one
gain_from_dist2 = @(d2) cosd(sqrt(d2)).^gain_exp;
gain_from_dist2_dB = @(d2) 10*log10(gain_from_dist2(d2));

angle_from_dist2_deg = @(d2) asind(min(1, sqrt(max(0, d2))));


main_beam_idx = nan(Cfg.NumUEs, nT);
main_beam_signal_lin = nan(Cfg.NumUEs, nT);
interference_lin = nan(Cfg.NumUEs, nT);
sir_lin = nan(Cfg.NumUEs, nT);
sir_dB = nan(Cfg.NumUEs, nT);

for i = 1:Cfg.NumUEs
    valid_t = find(served_mask(i,:));
    if isempty(valid_t)
        continue;
    end

    u = u_ues(i, valid_t)';
    v = v_ues(i, valid_t)';
    dist2_to_beams = (u - b_u').^2 + (v - b_v').^2;
    [min_dist2, mb_idx_local] = min(dist2_to_beams, [], 2);

    inside_main_beam = min_dist2 <= r_beam^2;
    if ~any(inside_main_beam)
        continue;
    end

    % Filter down to only times when the UE is actually inside a beam
    valid_t = valid_t(inside_main_beam);
    u = u(inside_main_beam);
    v = v(inside_main_beam);
    min_dist2 = min_dist2(inside_main_beam);
    mb_idx_local = mb_idx_local(inside_main_beam);

    main_beam_idx(i, valid_t) = mb_idx_local;
    
    % 1. Calculate Main Beam Signal
    sig_lin = gain_from_dist2(min_dist2);

    % ===============================================================
    % THE VECTORIZED NEIGHBOR CALCULATION
    % ===============================================================
    
    % Get an N_time x 6 matrix of neighbor indices for our serving beams
    nbs_matrix = neighbor_idx(mb_idx_local, :); 
    
    % Create a logical mask of where neighbors actually exist (not NaN)
    valid_nbs_mask = ~isnan(nbs_matrix);
    
    % Temporarily replace NaNs with 1s so we can use them as matrix indices without crashing
    safe_nbs_matrix = nbs_matrix;
    safe_nbs_matrix(~valid_nbs_mask) = 1;
    
    % Extract the u and v coordinates of all 6 neighbors simultaneously
    % These become N_time x 6 matrices
    b_u_nbs = b_u(safe_nbs_matrix); 
    b_v_nbs = b_v(safe_nbs_matrix);
    
    % Calculate distance squared to all 6 neighbors at once
    % u and v are N_time x 1, but MATLAB implicitly expands them to N_time x 6
    d2_nb = (u - b_u_nbs).^2 + (v - b_v_nbs).^2;
    
    % Run the gain function on the entire N_time x 6 matrix
    int_gains = gain_from_dist2(d2_nb);
    
    % Zero out the gains for the "fake" neighbors we mapped to index 1 earlier
    int_gains(~valid_nbs_mask) = 0;
    
    % Sum horizontally across the 6 columns to get total interference per timestep
    int_lin = sum(int_gains, 2);
    
    % ===============================================================

    % Calculate final SIR
    sir_val_lin = sig_lin ./ max(int_lin, eps);

    main_beam_signal_lin(i, valid_t) = sig_lin;
    interference_lin(i, valid_t) = int_lin;
    sir_lin(i, valid_t) = sir_val_lin;
    sir_dB(i, valid_t) = 10*log10(sir_val_lin);
end

Interference = struct();
Interference.main_beam_idx = main_beam_idx;
Interference.main_beam_signal_lin = main_beam_signal_lin;
Interference.interference_lin = interference_lin;
Interference.SIR_lin = sir_lin;
Interference.SIR_dB = sir_dB;
Interference.valid_samples = isfinite(sir_dB);

Interference.mean_SIR_dB = mean(sir_dB, 2, 'omitnan');
Interference.p10_SIR_dB = nan(Cfg.NumUEs, 1);
for i = 1:Cfg.NumUEs
    this_sir = sir_dB(i, :);
    this_sir = this_sir(isfinite(this_sir));
    if ~isempty(this_sir)
        Interference.p10_SIR_dB(i) = prctile(this_sir, 10);
    end
end

all_sir = sir_dB(isfinite(sir_dB));
if isempty(all_sir)
    Interference.global_mean_SIR_dB = NaN;
    Interference.global_p10_SIR_dB = NaN;
else
    Interference.global_mean_SIR_dB = mean(all_sir);
    Interference.global_p10_SIR_dB = prctile(all_sir, 10);
end

fprintf('Interference model complete. Global mean SIR = %.2f dB, p10 SIR = %.2f dB\n', ...
    Interference.global_mean_SIR_dB, Interference.global_p10_SIR_dB);

%% 6. Diagnostic Plots to Validate Interference Calculation
% Pick one UE/time sample where the serving beam has all 6 neighbors
candidate_u = [];
candidate_t = [];
candidate_main_dist2 = [];

for i = 1:Cfg.NumUEs
    valid_t = find(isfinite(main_beam_idx(i,:)));
    if isempty(valid_t)
        continue;
    end

    for tt = valid_t
        mb = main_beam_idx(i, tt);
        if ~isfinite(mb)
            continue;
        end
        nbs = neighbor_idx(mb, :);
        if sum(isfinite(nbs)) ~= 6
            continue;
        end

        d2_main = (u_ues(i,tt) - b_u(mb)).^2 + (v_ues(i,tt) - b_v(mb)).^2;
        candidate_u(end+1,1) = i; %#ok<SAGROW>
        candidate_t(end+1,1) = tt; %#ok<SAGROW>
        candidate_main_dist2(end+1,1) = d2_main; %#ok<SAGROW>
    end
end

if isempty(candidate_u)
    % Fallback: any valid sample
    [u_idx, t_idx] = find(isfinite(main_beam_idx), 1, 'first');
else
    [~, best_idx] = min(candidate_main_dist2);
    u_idx = candidate_u(best_idx);
    t_idx = candidate_t(best_idx);
end

if ~isempty(u_idx)
    mb = main_beam_idx(u_idx, t_idx);
    nbs = neighbor_idx(mb, :);
    nbs = nbs(isfinite(nbs));
    beam_set = [mb, nbs];

    u0 = u_ues(u_idx, t_idx);
    v0 = v_ues(u_idx, t_idx);
    d2_set = (u0 - b_u(beam_set)).^2 + (v0 - b_v(beam_set)).^2;
    gain_set_lin = gain_from_dist2(d2_set);
    gain_set_dB = gain_from_dist2_dB(d2_set);
    theta_set_deg = angle_from_dist2_deg(d2_set);

    sig_snap = gain_set_lin(1);
    int_snap = sum(gain_set_lin(2:end));
    sir_snap_dB = 10*log10(sig_snap / max(int_snap, eps));
    sir_stored_dB = sir_dB(u_idx, t_idx);

    out_dir = fullfile('screenshots', 'interference_debug');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % Plot A: Beam geometry + gains at one instant
    f_diag1 = figure('Color', 'w', 'Position', [80 80 1300 520]);
    tlo = tiledlayout(f_diag1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(tlo, 1); hold on; grid on; box on; axis equal;
    theta = linspace(0, 2*pi, 200);

    for k = 1:numel(beam_set)
        bi = beam_set(k);
        x = b_u(bi) + r_beam*cos(theta);
        y = b_v(bi) + r_beam*sin(theta);
        if k == 1
            plot(x, y, 'b-', 'LineWidth', 2.2);
            scatter(b_u(bi), b_v(bi), 70, 'b', 'filled');
        else
            plot(x, y, '-', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);
            scatter(b_u(bi), b_v(bi), 50, [0.85 0.2 0.2], 'filled');
        end

        txt = sprintf('G=%.2f dB | %.2f deg', gain_set_dB(k), theta_set_deg(k));
        text(b_u(bi), b_v(bi) + 0.01, txt, 'HorizontalAlignment', 'center', ...
            'FontSize', 9, 'FontWeight', 'bold');
        plot([u0 b_u(bi)], [v0 b_v(bi)], ':', 'Color', [0.45 0.45 0.45]);
    end

    scatter(u0, v0, 120, 'k', 'filled');
    text(u0, v0 - 0.012, sprintf('UE %d', u_idx), 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 10);

    title(sprintf('7-Beam Snapshot (UE %d, t=%d)', u_idx, t_idx), 'FontWeight', 'bold');
    xlabel('u = sin(\eta)cos(\phi)');
    ylabel('v = sin(\eta)sin(\phi)');
    legend({'Main beam','Main center','Neighbor beam','Neighbor center','Link to UE','UE'}, ...
        'Location', 'southoutside');

    nexttile(tlo, 2); hold on; grid on; box on;
    labels = strings(1, numel(beam_set));
    labels(1) = "Main";
    for k = 2:numel(beam_set)
        labels(k) = "N" + string(k-1);
    end
    bar(categorical(labels), gain_set_dB, 'FaceColor', [0.2 0.5 0.85]);
    ylabel('Gain (dB, relative)');
    title(sprintf('Main/Neighbor Gains | SIR_{snap}=%.2f dB (stored %.2f dB)', ...
        sir_snap_dB, sir_stored_dB), 'FontWeight', 'bold');

    exportgraphics(f_diag1, fullfile(out_dir, 'Beam_Gain_Snapshot.png'), 'Resolution', 300);

    % Plot B: SIR across time for the selected UE
    f_diag2 = figure('Color', 'w', 'Position', [120 120 1200 420]);
    sir_series = sir_dB(u_idx, :);
    % time_vec = metrics.SimData(u_idx).Time;
    time_vec = minutes(metrics.SimData(u_idx).Time - metrics.SimData(u_idx).Time(1));

    if numel(time_vec) == nT
        plot(time_vec, sir_series, 'LineWidth', 1.3, 'Color', [0 0.45 0.74]);
        xlabel('Time');
    else
        plot(1:nT, sir_series, 'LineWidth', 1.3, 'Color', [0 0.45 0.74]);
        xlabel('Time index');
    end

    grid on; box on;
    ylabel('SIR (dB)');
    title(sprintf('SIR Over Time for UE %d | Mean=%.2f dB, P10=%.2f dB', ...
        u_idx, Interference.mean_SIR_dB(u_idx), Interference.p10_SIR_dB(u_idx)), ...
        'FontWeight', 'bold');

    exportgraphics(f_diag2, fullfile(out_dir, 'SIR_Time_Series_UE.png'), 'Resolution', 300);

    % Plot C: Main signal, neighbor interference, and SIR consistency over time
    f_diag3 = figure('Color', 'w', 'Position', [130 130 1300 650]);
    t3 = tiledlayout(f_diag3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    sig_lin_series = main_beam_signal_lin(u_idx, :);
    int_lin_series = interference_lin(u_idx, :);
    sir_recomputed_dB = 10*log10(sig_lin_series ./ max(int_lin_series, eps));
    valid_series = isfinite(sir_dB(u_idx, :)) & isfinite(sir_recomputed_dB);

    sig_dB_series = 10*log10(max(sig_lin_series, eps));
    int_dB_series = 10*log10(max(int_lin_series, eps));

    if numel(time_vec) == nT
        xvals = time_vec;
        xlab = 'Time';
    else
        xvals = 1:nT;
        xlab = 'Time index';
    end

    nexttile(t3, 1); hold on; grid on; box on;
    plot(xvals, sig_dB_series, 'LineWidth', 1.2, 'Color', [0 0.45 0.74]);
    plot(xvals, int_dB_series, 'LineWidth', 1.2, 'Color', [0.85 0.33 0.1]);
    ylabel('Relative Gain (dB)');
    xlabel(xlab);
    title(sprintf('UE %d: Main Signal vs Neighbor Interference', u_idx), 'FontWeight', 'bold');
    legend({'Main signal', 'Sum neighbor interference'}, 'Location', 'best');

    nexttile(t3, 2); hold on; grid on; box on;
    plot(xvals, sir_dB(u_idx, :), 'LineWidth', 1.3, 'Color', [0.49 0.18 0.56]);
    plot(xvals, sir_recomputed_dB, '--', 'LineWidth', 1.1, 'Color', [0.1 0.1 0.1]);
    ylabel('SIR (dB)');
    xlabel(xlab);

    if any(valid_series)
        max_abs_err = max(abs(sir_dB(u_idx, valid_series) - sir_recomputed_dB(valid_series)));
    else
        max_abs_err = NaN;
    end

    title(sprintf('SIR Consistency Check (max |Delta| = %.3e dB)', max_abs_err), 'FontWeight', 'bold');
    legend({'Stored SIR', 'Recomputed SIR'}, 'Location', 'best');

    exportgraphics(f_diag3, fullfile(out_dir, 'Signal_Interference_SIR_Consistency.png'), 'Resolution', 300);

    fprintf('Saved diagnostic plots to: %s\n', out_dir);
    fprintf('Snapshot consistency check: recomputed SIR = %.3f dB, stored SIR = %.3f dB\n', ...
        sir_snap_dB, sir_stored_dB);
    fprintf('Time-series consistency check for UE %d: max abs error = %.3e dB\n', u_idx, max_abs_err);
else
    warning('No valid UE/time sample found for diagnostic plotting.');
end