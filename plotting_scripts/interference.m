clearvars; close all; clc;

%% 1) Single-UE scenario (Aalborg)
constellation_type = "walkerdelta";
ue_grid_size = "small";
duration = "short";
frequency = "ku";
height_km = 1000;

Cfg = get_cfg(height_km, constellation_type, ue_grid_size, duration, frequency);
Cfg.SampleTime = 0.1; % seconds
Cfg.StopTime = datetime('1-Jun-2025 12:29:59', 'TimeZone', 'UTC');

% Force exactly one UE in Aalborg.
Cfg.Flat_UE_array.Lats = 57.0488;
Cfg.Flat_UE_array.Lons = 9.9217;

calc_link = true;
use_parallel = false;
metrics = coverage_simulator_function(Cfg, use_parallel, calc_link);

UE = metrics.UEs(1);
t = UE.SimData.Time;
el = UE.SimData.Elevation_deg;
az = UE.SimData.Azimuth_deg;
num_visible = UE.SimData.Num_visible;
sat_id = UE.SimData.SatID;

sir_dB = UE.DL.SIR;
sinr_dB = UE.DL.SINR;
snr_dB = UE.DL.SNR;
throughput_Mbps = UE.DL.Throughput / 1e6;

valid_sir = isfinite(sir_dB);
valid_service = isfinite(el);
if ~any(valid_service)
    warning('UE has no valid service/elevation samples in the selected window. Extend StopTime or lower Min_elevation_UE.');
    return;
end

% Choose best available quality metric for snapshot selection.
if any(valid_sir)
    metric_name = 'SIR';
    metric_series_dB = sir_dB;
elseif any(isfinite(sinr_dB))
    metric_name = 'SINR';
    metric_series_dB = sinr_dB;
else
    metric_name = 'SNR';
    metric_series_dB = snr_dB;
end

valid_metric = isfinite(metric_series_dB);
metric_valid = metric_series_dB(valid_metric);
mean_metric = mean(metric_valid);
p10_metric = prctile(metric_valid, 10);

[best_metric, best_rel_idx] = max(metric_valid);
[worst_metric, worst_rel_idx] = min(metric_valid);
valid_idx = find(valid_metric);
t_best = valid_idx(best_rel_idx);
t_worst = valid_idx(worst_rel_idx);

fprintf('Single UE (Aalborg) complete. Mean %s = %.2f dB, P10 %s = %.2f dB\n', metric_name, mean_metric, metric_name, p10_metric);
fprintf('Best %s = %.2f dB at %s\n', metric_name, best_metric, datestr(t(t_best)));
fprintf('Worst %s = %.2f dB at %s\n', metric_name, worst_metric, datestr(t(t_worst)));

%% 2) Plots from online link-calculation output
out_dir = fullfile('figures', 'interference');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end



% % B) Throughput and elevation
% f2 = figure('Color', 'w', 'Position', [130 130 1200 420]);
% yyaxis left;
% plot(t, throughput_Mbps, 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74]);
% ylabel('Throughput (Mbps)');
% yyaxis right;
% plot(t, el, 'LineWidth', 1.3, 'Color', [0.47 0.67 0.19]);
% ylabel('Elevation (deg)');
% grid on; box on;
% xlabel('Time');
% title('Aalborg UE Throughput and Elevation', 'FontWeight', 'bold');
% legend({'Throughput','Elevation'}, 'Location', 'best');
% exportgraphics(f2, fullfile(out_dir, 'Aalborg_Throughput_Elevation.png'), 'Resolution', 600);

% % C) Visibility and serving satellite ID
% f3 = figure('Color', 'w', 'Position', [140 140 1200 420]);
% yyaxis left;
% stairs(t, num_visible, 'LineWidth', 1.3, 'Color', [0.30 0.30 0.30]);
% ylabel('Visible satellites');
% yyaxis right;
% stairs(t, sat_id, 'LineWidth', 1.1, 'Color', [0.85 0.33 0.10]);
% ylabel('Serving Sat ID');
% grid on; box on;
% xlabel('Time');
% title(sprintf('Aalborg UE Visibility / Serving Satellite | Azimuth valid samples: %d', sum(isfinite(az))), 'FontWeight', 'bold');
% legend({'Num visible','Serving SatID'}, 'Location', 'best');
% exportgraphics(f3, fullfile(out_dir, 'Aalborg_Visibility_SatID.png'), 'Resolution', 600);

% fprintf('Saved current diagnostics to: %s\n', out_dir);

%% 3) Restore the old best/worst beam diagnostics
out_dir = fullfile('figures', 'interference');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% A) SIR/SINR/SNR time series
f1 = figure('Color', 'w', 'Position', [120 120 1200 420]);
hold on; grid on; box on;
plot(t, sir_dB, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
% plot(t, sinr_dB, 'LineWidth', 1.3, 'Color', [0.00 0.45 0.74]);
% plot(t, snr_dB, 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
scatter(t(t_best), best_metric, 70, 'g', 'filled', 'MarkerEdgeColor', 'k');
scatter(t(t_worst), worst_metric, 70, 'r', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Time');
ylabel('dB');
title(sprintf('Aalborg UE Link Quality | Mean %s = %.2f dB, P10 = %.2f dB', metric_name, mean_metric, p10_metric), 'FontWeight', 'bold');
legend({'SIR',['Best ' metric_name],['Worst ' metric_name]}, 'Location', 'best');
exportgraphics(f1, fullfile(out_dir, 'Aalborg_SIR_SINR_TimeSeries.png'), 'Resolution', 600);

% Reconstruct the beam grid used by the online link calculation.
BeamGrid = calculate_Beams(Cfg.DL.f, Cfg.DL.G_tx, Cfg.Orbit_height, Cfg.Min_elevation_UE, Cfg.FRF);
b_u = BeamGrid.b_u;
b_v = BeamGrid.b_v;
neighbor_idx = BeamGrid.neighbor_idx;
r_beam = BeamGrid.r_beam;

serving_beam_idx = UE.DL.serving_beam_idx;
serving_beam_signal_lin = UE.DL.serving_beam_signal_lin;
interference_lin = UE.DL.interference_lin;

sir_series = UE.DL.SIR;
sir_valid = isfinite(sir_series) & isfinite(serving_beam_idx) & serving_beam_idx > 0;

if ~any(sir_valid)
    warning('No valid SIR samples found for old-style diagnostics. Skipping restored plots.');
    return;
end

sir_valid_vals = sir_series(sir_valid);
mean_sir = mean(sir_valid_vals);
p10_sir = prctile(sir_valid_vals, 10);
[max_sir, best_rel_idx] = max(sir_valid_vals);
[min_sir, worst_rel_idx] = min(sir_valid_vals);
valid_idx = find(sir_valid);
t_best = valid_idx(best_rel_idx);
t_worst = valid_idx(worst_rel_idx);

time_cases = [t_best, t_worst];
case_names = {"Best Case", "Worst Case"};
file_names = {"Beam_Snapshot_BEST.png", "Beam_Snapshot_WORST.png"};

N_elements_snap = round(101.2 / BeamGrid.Beamwidth_deg);
d_spacing_snap = BeamGrid.lambda / 2;
phase_term_snap = @(d2) (pi * d_spacing_snap / BeamGrid.lambda) * max(sqrt(max(0, d2)), eps);
array_factor_snap = @(d2) sin(N_elements_snap * phase_term_snap(d2)) ./ (N_elements_snap * sin(phase_term_snap(d2)));
gain_from_dist2_snap = @(d2) array_factor_snap(d2).^2;
gain_from_dist2_dB_snap = @(d2) 10*log10(gain_from_dist2_snap(d2));

for c = 1:2
    t_idx = time_cases(c);

    mb = serving_beam_idx(1, t_idx);
    nbs = neighbor_idx(mb, :);
    nbs = nbs(isfinite(nbs));
    beam_set = [mb, nbs];

    u0 = sind(asind((6378.14e3 / (6378.14e3 + Cfg.Orbit_height)) * cosd(el(t_idx)))) * cosd(az(t_idx) + 180);
    v0 = sind(asind((6378.14e3 / (6378.14e3 + Cfg.Orbit_height)) * cosd(el(t_idx)))) * sind(az(t_idx) + 180);

    d2_set = (u0 - b_u(beam_set)).^2 + (v0 - b_v(beam_set)).^2;

    gain_set_lin = gain_from_dist2_snap(d2_set);
    gain_set_dB = gain_from_dist2_dB_snap(d2_set);

    sig_snap = gain_set_lin(1);
    gain_set_lin(2:end) = gain_set_lin(2:end) * Cfg.RU;
    gain_set_dB(2:end) = 10*log10(gain_set_lin(2:end));
    int_snap = sum(gain_set_lin(2:end));
    sir_snap_dB = 10*log10(sig_snap / max(int_snap, eps));
    sir_live_dB = sir_series(t_idx);
    sir_delta_dB = sir_snap_dB - sir_live_dB;

    f_diag1 = figure('Color', 'w', 'Position', [100 100 550 550]);
    hold on; grid on; box on; axis equal;
    theta = linspace(0, 2*pi, 200);
    h_serving_beam = [];
    h_serving_center = [];
    h_neighbor_beam = [];
    h_neighbor_center = [];
    h_link_to_ue = [];
    for k = 1:numel(beam_set)
        bi = beam_set(k);
        x = b_u(bi) + r_beam*cos(theta);
        y = b_v(bi) + r_beam*sin(theta);
        if k == 1
            h_serving_beam = plot(x, y, 'b-', 'LineWidth', 2.2);
            h_serving_center = scatter(b_u(bi), b_v(bi), 70, 'b', 'filled');
        else
            h_nb = plot(x, y, '-', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);
            h_nc = scatter(b_u(bi), b_v(bi), 50, [0.85 0.2 0.2], 'filled');
            if isempty(h_neighbor_beam)
                h_neighbor_beam = h_nb;
            end
            if isempty(h_neighbor_center)
                h_neighbor_center = h_nc;
            end
        end
        txt = sprintf('G=%.2f dB', gain_set_dB(k));
        text(b_u(bi), b_v(bi) + 0.01, txt, 'HorizontalAlignment', 'center', ...
            'FontSize', 9, 'FontWeight', 'bold');
        h_link = plot([u0 b_u(bi)], [v0 b_v(bi)], ':', 'Color', [0.45 0.45 0.45]);
        if isempty(h_link_to_ue)
            h_link_to_ue = h_link;
        end
    end
    h_ue = scatter(u0, v0, 120, 'k', 'filled');
    % text(u0, v0 - 0.012, sprintf('UE %d', 1), 'HorizontalAlignment', 'center', ...
    %     'FontWeight', 'bold', 'FontSize', 10);
    title(sprintf('%s (t=%d, RU=%.1f, FRF=%d)', case_names{c}, t_idx, Cfg.RU, Cfg.FRF), 'FontWeight', 'bold');
    xlabel('u = sin(\eta)cos(\phi)');
    ylabel('v = sin(\eta)sin(\phi)');
    legend_handles = [h_serving_beam, h_serving_center];
    legend_labels = {'Serving HPBW','Serving center'};
    if ~isempty(h_neighbor_beam)
        legend_handles = [legend_handles, h_neighbor_beam];
        legend_labels = [legend_labels, {'Neighbor HPBW'}];
    end
    if ~isempty(h_neighbor_center)
        legend_handles = [legend_handles, h_neighbor_center];
        legend_labels = [legend_labels, {'Neighbor center'}];
    end
    legend_handles = [legend_handles, h_link_to_ue, h_ue];
    legend_labels = [legend_labels, {'Link to UE','UE'}];
    legend(legend_handles, legend_labels, 'Location', 'southoutside', 'NumColumns', 3);
    % text(0.02, 0.98, sprintf('SIR = %.2f dB', sir_snap_dB), 'Units', 'normalized', ...
        % 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold');

    exportgraphics(f_diag1, fullfile(out_dir, file_names{c}), 'Resolution', 600);
    close(f_diag1);
end

time_vec = UE.SimData.Time;

f_diag2 = figure('Color', 'w', 'Position', [100 100 1000 550]);
hold on; grid on; box on;
plot(time_vec, sir_series, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
scatter(time_vec(t_best), max_sir, 100, 'g', 'filled', 'MarkerEdgeColor', 'k');
scatter(time_vec(t_worst), min_sir, 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Time');
ylabel('SIR (dB)');
title(sprintf('SIR Over Time for UE %d | Mean = %.2f dB, P10 = %.2f dB, RU=%.1f, FRF=%d', ...
    1, mean_sir, p10_sir, Cfg.RU, Cfg.FRF), 'FontWeight', 'bold');
legend({'SIR', 'Best Case Snapshot', 'Worst Case Snapshot'}, 'Location', 'best');
exportgraphics(f_diag2, fullfile(out_dir, 'SIR_Time_Series_UE.png'), 'Resolution', 600);

sig_lin_series = serving_beam_signal_lin(1, :);
int_lin_series = interference_lin(1, :);
sig_dB_series = 10*log10(max(sig_lin_series, eps));
int_dB_series = 10*log10(max(int_lin_series, eps));

f_diag3 = figure('Color', 'w', 'Position', [100 100 1000 550]);
hold on; grid on; box on;
plot(time_vec, sig_dB_series, 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
plot(time_vec, int_dB_series, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]);
xlabel('Time');
ylabel('Relative Gain (dB)');
end_idx = min(1200, numel(time_vec));
xlim([time_vec(1), time_vec(end_idx)]);
title(sprintf('UE %d: Serving Signal vs. Sum Neighbor Interference, RU=%.1f, FRF=%d', 1, Cfg.RU, Cfg.FRF), 'FontWeight', 'bold');
legend({'Serving Signal', 'Total Interference'}, 'Location', 'best');
exportgraphics(f_diag3, fullfile(out_dir, 'Signal_vs_Interference.png'), 'Resolution', 600);

f = Cfg.DL.f;
G = Cfg.DL.G_tx;
Beamwidth_deg = BeamGrid.Beamwidth_deg;
N_elements = round(101.2 / Beamwidth_deg);
% With d = lambda/2 spacing, the normalized array-factor shape is frequency-invariant.
phase_term = @(d2) (pi/2) * max(sqrt(max(0, d2)), eps);
array_factor = @(d2) sin(N_elements * phase_term(d2)) ./ (N_elements * sin(phase_term(d2)));
gain_from_dist2_dB = @(d2) 10*log10(array_factor(d2).^2);

theta_deg = linspace(-10, 10, 2000);
d2_1D = sind(theta_deg).^2;
beam_pattern_dB = gain_from_dist2_dB(d2_1D);

du = sind(Beamwidth_deg) / 2;
angle_frf1 = asind(sqrt(3) * du);
angle_frf3 = asind(3 * du);

figure('Name', 'Phased Array Beam Pattern', 'Color', 'w', 'Position', [100 100 1000 550]);
hold on; grid on; box on;
plot(theta_deg, beam_pattern_dB, 'LineWidth', 2, 'Color', '#0072BD');
yline(-3, '--', '-3 dB HPBW Crossover', 'Color', '#777777', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
xline(angle_frf1, '-.', 'FRF=1 Neighbor beam center', 'Color', '#D95319', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(-angle_frf1, '-.', 'Color', '#D95319', 'LineWidth', 1.5);
xline(angle_frf3, '-.', 'FRF=3 Neighbor beam center', 'Color', '#EDB120', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(-angle_frf3, '-.', 'Color', '#EDB120', 'LineWidth', 1.5);
ylim([-45 5]);
xlim([-10 10]);
xlabel('Steering Angle \theta (deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Gain [dB]', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Phased Array Radiation Pattern (%d Elements)\nServing Lobe vs. Neighboring Beam Locations', N_elements), 'FontSize', 14);
exportgraphics(gcf, fullfile(out_dir, 'Antenna_Pattern_1D.png'), 'Resolution', 600);

fprintf(sprintf('SIR = %.2f dB (live %.2f, dSIR %.3f dB)', sir_snap_dB, sir_live_dB, sir_delta_dB));

fprintf('Saved plots: %s\n', out_dir);