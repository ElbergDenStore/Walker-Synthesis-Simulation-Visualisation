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
out_dir = fullfile('plotting_scripts/figures', 'interference');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% A) SIR/SINR/SNR time series
f1 = figure('Color', 'w', 'Position', [120 120 500 250]);
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
exportgraphics(f1, fullfile(out_dir, 'Aalborg_SIR_SINR_TimeSeries.png'), 'Resolution', 300);

% Reconstruct the beam grid used by the online link calculation.
if isfield(Cfg, 'DL') && isfield(Cfg.DL, 'BeamGrid') && ~isempty(Cfg.DL.BeamGrid)
    BeamGrid = Cfg.DL.BeamGrid;
else
    error('Plotting requires Cfg.DL.BeamGrid to be populated before the simulation is run.');
end
b_u = BeamGrid.u_center;
b_v = BeamGrid.v_center;
neighbor_idx = BeamGrid.neighbor_idx;
r_beam = 1.391 / (BeamGrid.Nu * (pi/2));

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

% Gain function using the 2D separable array factor (matches beam_gain_and_interference.m)
Nu_snap = BeamGrid.Nu;
Nv_snap = BeamGrid.Nv;

% Pre-compute the exact AF-only -3 dB isoline template (radial sweep + binary search).
% The true 2D -3 dB contour of the separable pattern is a rounded square, NOT a circle.
% This template is centred at origin and gets translated to each beam centre in the loop.
theta_iso = linspace(0, 2*pi, 360);
r_iso = zeros(1, 360);
for ki = 1:360
    ct = cos(theta_iso(ki)); st = sin(theta_iso(ki));
    r_lo = 0; r_hi = 4 * r_beam;
    for iter = 1:40
        rm = (r_lo + r_hi) / 2;
        du_t = rm * ct; if du_t == 0, du_t = eps; end
        dv_t = rm * st; if dv_t == 0, dv_t = eps; end
        af_u = sin(Nu_snap*(pi/2)*du_t) / (Nu_snap*sin((pi/2)*du_t));
        af_v = sin(Nv_snap*(pi/2)*dv_t) / (Nv_snap*sin((pi/2)*dv_t));
        if 20*log10(abs(af_u*af_v)) > -3
            r_lo = rm;
        else
            r_hi = rm;
        end
    end
    r_iso(ki) = (r_lo + r_hi) / 2;
end
iso_x_tmpl = r_iso .* cos(theta_iso);
iso_y_tmpl = r_iso .* sin(theta_iso);

for c = 1:2
    t_idx = time_cases(c);

    mb = serving_beam_idx(1, t_idx);
    nbs = neighbor_idx(mb, :);
    nbs = nbs(isfinite(nbs));
    % First ring only: all nearest co-channel beams share the same distance in a hex grid.
    if ~isempty(nbs)
        d_nbs = sqrt((b_u(nbs) - b_u(mb)).^2 + (b_v(nbs) - b_v(mb)).^2);
        nbs = nbs(d_nbs <= d_nbs(1) * 1.05);
    end
    beam_set = [mb, nbs];

    u0 = sind(asind((6378.14e3 / (6378.14e3 + Cfg.Orbit_height)) * cosd(el(t_idx)))) * cosd(az(t_idx) + 180);
    v0 = sind(asind((6378.14e3 / (6378.14e3 + Cfg.Orbit_height)) * cosd(el(t_idx)))) * sind(az(t_idx) + 180);

    du_set = u0 - b_u(beam_set);
    dv_set = v0 - b_v(beam_set);
    du_set(du_set == 0) = eps;
    dv_set(dv_set == 0) = eps;
    AF_u_set = sin(Nu_snap * (pi/2) .* du_set) ./ (Nu_snap .* sin((pi/2) .* du_set));
    AF_v_set = sin(Nv_snap * (pi/2) .* dv_set) ./ (Nv_snap .* sin((pi/2) .* dv_set));
    gain_set_dB = 20 * log10(abs(AF_u_set .* AF_v_set));

    % Add relative element factor (flat phased array: nadir-pointing face).
    % Matches calculate_total_gain_dB in beam_gain_and_interference.m.
    EF_ue_dB       = 10 * BeamGrid.Cos_exponent * log10(max(sqrt(1 - min(u0^2 + v0^2, 1)), eps));
    rho_sq_centers = min(b_u(beam_set).^2 + b_v(beam_set).^2, 1);
    EF_centers_dB  = 10 * BeamGrid.Cos_exponent * log10(max(sqrt(1 - rho_sq_centers), eps));
    gain_set_dB    = gain_set_dB + (EF_ue_dB - EF_centers_dB);

    gain_set_lin = 10.^(gain_set_dB / 10);

    sig_snap = gain_set_lin(1);
    gain_set_lin(2:end) = gain_set_lin(2:end) * Cfg.RU;
    gain_set_dB(2:end) = 10*log10(gain_set_lin(2:end));
    int_snap = sum(gain_set_lin(2:end));
    sir_snap_dB = 10*log10(sig_snap / max(int_snap, eps));
    sir_live_dB = sir_series(t_idx);
    sir_delta_dB = sir_snap_dB - sir_live_dB;

    f_diag1 = figure('Color', 'w', 'Position', [100 100 600 400]);
    hold on; grid on; box on; axis equal;
    h_serving_beam = [];
    h_serving_center = [];
    h_neighbor_beam = [];
    h_neighbor_center = [];
    h_link_to_ue = [];
    for k = 1:numel(beam_set)
        bi = beam_set(k);
        x = b_u(bi) + iso_x_tmpl;
        y = b_v(bi) + iso_y_tmpl;
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

    exportgraphics(f_diag1, fullfile(out_dir, file_names{c}), 'Resolution', 300);
    close(f_diag1);
end

time_vec = UE.SimData.Time;

f_diag2 = figure('Color', 'w', 'Position', [100 100 500 500]);
hold on; grid on; box on;
plot(time_vec, sir_series, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
scatter(time_vec(t_best), max_sir, 100, 'g', 'filled', 'MarkerEdgeColor', 'k');
scatter(time_vec(t_worst), min_sir, 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Time');
ylabel('SIR (dB)');
title(sprintf('SIR Over Time for UE %d | Mean = %.2f dB, P10 = %.2f dB, RU=%.1f, FRF=%d', ...
    1, mean_sir, p10_sir, Cfg.RU, Cfg.FRF), 'FontWeight', 'bold');
legend({'SIR', 'Best Case Snapshot', 'Worst Case Snapshot'}, 'Location', 'best');
exportgraphics(f_diag2, fullfile(out_dir, 'SIR_Time_Series_UE.png'), 'Resolution', 300);

% serving_beam_signal_lin stores EIRP density in dBm/Hz (already in dB)
% interference_lin stores sum interference power density in linear mW/Hz
sig_dB_series = serving_beam_signal_lin(1, :);
int_dB_series  = 10*log10(max(interference_lin(1, :), eps));

f_diag3 = figure('Color', 'w', 'Position', [100 100 400 400]);
hold on; grid on; box on;
plot(time_vec, sig_dB_series, 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
plot(time_vec, int_dB_series, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]);
xlabel('Time');
ylabel('Relative Gain (dB)');
end_idx = min(1200, numel(time_vec));
xlim([time_vec(1), time_vec(end_idx)]);
title(sprintf('UE %d: Serving Signal vs. Sum Neighbor Interference, RU=%.1f, FRF=%d', 1, Cfg.RU, Cfg.FRF), 'FontSize', 14);
legend({'Serving Signal', 'Total Interference'}, 'Location', 'best');
exportgraphics(f_diag3, fullfile(out_dir, 'Signal_vs_Interference.png'), 'Resolution', 300);

f = Cfg.DL.f;
G = Cfg.DL.G_tx;
N_elements = BeamGrid.Nu;
% With d = lambda/2 spacing, the normalized array-factor shape is frequency-invariant.
phase_term = @(d2) (pi/2) * max(sqrt(max(0, d2)), eps);
array_factor = @(d2) sin(N_elements * phase_term(d2)) ./ (N_elements * sin(phase_term(d2)));
gain_from_dist2_dB = @(d2) 10*log10(array_factor(d2).^2);

theta_deg = linspace(-10, 10, 2000);
d2_1D = sind(theta_deg).^2;
beam_pattern_dB = gain_from_dist2_dB(d2_1D);

% r_beam is the 3dB half-beamwidth radius in u/v space
angle_frf1 = asind(sqrt(3) * r_beam);
angle_frf3 = asind(3 * r_beam);

figure('Name', 'Phased Array Beam Pattern', 'Color', 'w', 'Position', [100 100 500 300]);
hold on; grid on; box on;
plot(theta_deg, beam_pattern_dB, 'LineWidth', 2, 'Color', '#0072BD');
yline(-3, '--', '-3 dB HPBW Crossover', 'Color', '#777777', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
xline(angle_frf1, '-.', 'FRF=1 Neighbor beam center', 'Color', '#D95319', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(-angle_frf1, '-.', 'Color', '#D95319', 'LineWidth', 1.5);
xline(angle_frf3, '-.', 'FRF=3 Neighbor beam center', 'Color', '#EDB120', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(-angle_frf3, '-.', 'Color', '#EDB120', 'LineWidth', 1.5);
ylim([-45 5]);
xlim([-10 10]);
xlabel('Steering Angle \theta (deg)', 'FontSize', 14);
ylabel('Relative Gain (dB)', 'FontSize', 14);
title(sprintf('Phased Array Radiation Pattern (%d Elements)\nServing Lobe vs. Neighboring Beam Locations', N_elements), 'FontSize', 14);
exportgraphics(gcf, fullfile(out_dir, 'Antenna_Pattern_1D.png'), 'Resolution', 300);

fprintf(sprintf('SIR = %.2f dB (live %.2f, dSIR %.3f dB)', sir_snap_dB, sir_live_dB, sir_delta_dB));

fprintf('Saved plots: %s\n', out_dir);