% PLOT_BEAM_UTILIZATION  Illustrate how user equipment spreads across the beams.
%   Runs a scenario and plots per-beam utilisation, showing which beams carry
%   load and how evenly traffic is distributed across the footprint.
clearvars; close all; clc;

%% Output Directory
out_dir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'plotting_scripts/figures', 'beam_utilization');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% Figure Style Configuration
% Edit these values to control all figure sizes, fonts, and export resolution.
FIG.dpi        = 300;       % export DPI
FIG.font_size  = 14;        % axis label / tick font size
FIG.title_size = 14;        % title font size (via multiplier)
FIG.annot_size = 9;         % small annotation / beam-label text
FIG.lw         = 1.5;       % default line width (where not explicitly overridden)
FIG.wide    = [100 100 600 400];
FIG.medium  = [100 100  500 400];
FIG.square  = [100 100  400 400];
FIG.compact = [100 100  300 300];
FIG.tall    = [100 100  400 400];
set(groot, 'DefaultAxesFontSize',                FIG.font_size);
set(groot, 'DefaultTextFontSize',                FIG.font_size);
set(groot, 'DefaultAxesTitleFontSizeMultiplier', FIG.title_size / FIG.font_size);
set(groot, 'DefaultLineLineWidth',               FIG.lw);

%% Configuration & Region Selection
% Toggle between 'Nordjylland', 'Denmark', 'Full3000', or 'Full30000'
REGION = 'Full3000';

switch REGION
    case 'Nordjylland'
        latlim = [56.5, 58.0];
        lonlim = [8.0, 11.0];
        people_per_ue = 300;
    case 'Denmark'
        latlim = [54.5, 58.0];
        lonlim = [8.0, 15.5];
        people_per_ue = 3000;
    case 'Full3000'
        latlim = [54+(35/60), 83+(40/60)];
        lonlim = [-(73+(10/60)), 33+(30/60)];
        people_per_ue = 3000;
    case 'Full30000'
        latlim = [54+(35/60), 83+(40/60)];
        lonlim = [-(73+(10/60)), 33+(30/60)];
        people_per_ue = 30000;
end

fprintf('Running simulation for %s...\n', REGION);
Cfg = default_config(1000,"walkerdelta");
calc_link = false; use_parallel = false;
[Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, Total_Pop] = generate_population_ues(latlim, lonlim, people_per_ue);
Cfg.NumUEs = length(Cfg.Flat_UE_array.Lats);
metrics = Constellation_simulator(Cfg, use_parallel, calc_link);

%% System-Wide Utilization Analysis
nT = length(metrics.SimData(1).Time);
TotalSats = Cfg.Total_sats;
time_vec = metrics.SimData(1).Time;
time_mins = minutes(time_vec - time_vec(1));

% Map which Sat each UE is using at each second
all_SatIDs = nan(Cfg.NumUEs, nT);
for i = 1:Cfg.NumUEs
    all_SatIDs(i, :) = metrics.SimData(i).SatID;
end

% Calculate SatUtilization [TotalSats x nT]
SatUtilization = zeros(TotalSats, nT);
for t = 1:nT
    active_sats = all_SatIDs(:, t);
    active_sats = active_sats(~isnan(active_sats));
    if ~isempty(active_sats)
        SatUtilization(:, t) = histcounts(active_sats, 0.5:(TotalSats+0.5));
    end
end

%% 3D Surface: Constellation Utilization vs Time
figure('Color', 'w', 'Name', 'Constellation Load', 'Position', FIG.wide);
[TimeGrid, SatGrid] = meshgrid(time_mins, 1:TotalSats);

s = surf(TimeGrid, SatGrid, SatUtilization);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
colormap(turbo);
view(-35, 60);
grid on;

cb = colorbar;
ylabel(cb, 'Connected UEs');
xlabel('Time (minutes from start)');
ylabel('Satellite ID');
zlabel('UE Count per Satellite');
% title(sprintf('%s: System-Wide Satellite Utilization', REGION));
title("System-Wide Satellite Utilization");
exportgraphics(gcf, fullfile(out_dir, 'constellation_surf.png'), 'Resolution', FIG.dpi);

%% Plot: Constellation Utilization (Active Satellites Only)
figure('Color', 'w', 'Name', 'Active Satellites Load', 'Position', FIG.wide);

active_sat_indices = find(max(SatUtilization, [], 2) > 0);
imagesc(time_mins, 1:TotalSats, SatUtilization);

xlabel('Time (minutes from start)');
ylabel('Satellite ID (Constellation Index)');
title("Utilization of Satellites");

cb = colorbar;
ylabel(cb, 'Number of Connected UEs');
set(gca, 'YDir', 'normal');
grid on;

exportgraphics(gcf, fullfile(out_dir, 'ActiveSatUtilization.png'), 'Resolution', FIG.dpi);


% --- Find the Peak Satellite and Time ---
[max_val, linear_idx] = max(SatUtilization(:));
[peakSat, t_peak] = ind2sub(size(SatUtilization), linear_idx);

fprintf('Peak Found: Sat %d at %s with %d UEs.\n', peakSat, datestr(time_vec(t_peak)), max_val);

%% Geometry & Beam Grid
Re = 6371; h = Cfg.Orbit_height/1000; Min_Elev_deg = 20;
f = 20e9; c = 3e8; lambda = c/f; G = 40;
Beamwidth_deg = sqrt(32400./(10.^(G/10)));
r_beam = sind(Beamwidth_deg / 2);
eta_max = asind((Re / (Re + h)) * cosd(Min_Elev_deg));
du = sind(Beamwidth_deg) / 2;
rings = ceil(sind(eta_max) / du) + 2;

b_u = []; b_v = [];
b_q = []; b_r = [];

for q = -rings:rings
    for r = -rings:rings
        u_val = du * sqrt(3) * (q + r/2);
        v_val = du * 1.5 * r;
        if (u_val^2 + v_val^2) <= sind(eta_max)^2
            b_u = [b_u; u_val];
            b_v = [b_v; v_val];
            b_q = [b_q; q];
            b_r = [b_r; r];
        end
    end
end

[~, sort_idx] = sortrows([b_v, b_u], 'descend');
b_u = b_u(sort_idx);
b_v = b_v(sort_idx);
b_q = b_q(sort_idx);
b_r = b_r(sort_idx);

num_beams = length(b_u);

%% Vectorize UE Projections (filtered by peak sat)
u_ues_peak = nan(Cfg.NumUEs, nT);
v_ues_peak = nan(Cfg.NumUEs, nT);

for i = 1:Cfg.NumUEs
    is_connected = (all_SatIDs(i, :) == peakSat);
    if any(is_connected)
        el = metrics.SimData(i).Elevation_deg(is_connected);
        az = metrics.SimData(i).Azimuth_deg(is_connected);
        eta_vec = asind((Re / (Re + h)) * cosd(el));
        u_ues_peak(i, is_connected) = sind(eta_vec) .* cosd(az + 180);
        v_ues_peak(i, is_connected) = sind(eta_vec) .* sind(az + 180);
    end
end

%% Resource Analysis for the Single Peak Satellite
beam_activity_matrix = false(num_beams, nT);
for t = 1:nT
    u_t = u_ues_peak(~isnan(u_ues_peak(:,t)), t);
    v_t = v_ues_peak(~isnan(v_ues_peak(:,t)), t);
    if ~isempty(u_t)
        for b = 1:num_beams
            if any((u_t - b_u(b)).^2 + (v_t - b_v(b)).^2 <= r_beam^2)
                beam_activity_matrix(b, t) = true;
            end
        end
    end
end

%% Snapshot at peak utilization (in degrees)
fig1 = figure('Color', 'w', 'Visible', 'off', 'Position', FIG.square); hold on; axis equal; grid on;
theta = linspace(0, 2*pi, 100);

b_eta = asind(sqrt(b_u.^2 + b_v.^2)); b_az = atan2d(b_v, b_u);
for b = 1:num_beams
    plot(b_eta(b)*cosd(b_az(b)) + (Beamwidth_deg/2)*cos(theta), ...
         b_eta(b)*sind(b_az(b)) + (Beamwidth_deg/2)*sin(theta), 'Color', [0.85 0.85 0.85]);
end

u_snap = u_ues_peak(:, t_peak); v_snap = v_ues_peak(:, t_peak);
ue_eta = asind(sqrt(u_snap.^2 + v_snap.^2)); ue_az = atan2d(v_snap, u_snap);
scatter(ue_eta.*cosd(ue_az), ue_eta.*sind(ue_az), 3, 'r');

plot(eta_max*cos(theta), eta_max*sin(theta), 'k--', 'LineWidth', 2);
title(sprintf('Sat %d Peak Load (%d/%d UEs)\nTime: %s', peakSat, max_val, Cfg.NumUEs, datestr(time_vec(t_peak))));
xlabel('Degrees off Nadir'); ylabel('Degrees off Nadir');
exportgraphics(fig1, fullfile(out_dir, 'peak_sat_fov.png'), 'Resolution', FIG.dpi);

fprintf('Done! Saved peak analysis for Sat %d.\n', peakSat);


%% Resource Load Analysis & Beam Assignment (Single Sat Perspective)
active_beams_time = zeros(nT, 1);
max_ues_per_beam_time = zeros(nT, 1);
beam_ue_count_matrix = zeros(num_beams, nT);
r_beam_sq = r_beam^2;

fprintf('Calculating resource load for Peak Sat %d...\n', peakSat);

for t = 1:nT
    u_t = u_ues_peak(~isnan(u_ues_peak(:,t)), t);
    v_t = v_ues_peak(~isnan(v_ues_peak(:,t)), t);

    if ~isempty(u_t)
        dist_sq_matrix = (u_t - b_u').^2 + (v_t - b_v').^2;
        [~, best_beam_idx] = min(dist_sq_matrix, [], 2);

        if ~isempty(best_beam_idx)
            temp_beam_counts = histcounts(best_beam_idx, 0.5:(num_beams+0.5));
            active_beams_time(t) = sum(temp_beam_counts > 0);
            max_ues_per_beam_time(t) = max(temp_beam_counts);
            beam_ue_count_matrix(:, t) = temp_beam_counts';
        end
    end
end

%% Plotting: The Resource Load (Active Beams)
figA = figure('Color', 'w', 'Position', FIG.medium);
tiledlayout(2,1, 'TileSpacing', 'compact');

nexttile;
plot(time_mins, active_beams_time, 'LineWidth', 2, 'Color', [0 0.447 0.741]);
grid on; ylabel('Active Beams');
title(sprintf('Sat %d Resource Profile: Active Beams', peakSat));
yline(mean(active_beams_time(active_beams_time>0)), '--k', 'Avg Active');

nexttile;
plot(time_mins, max_ues_per_beam_time, 'LineWidth', 2, 'Color', [0.85 0.32 0.1]);
grid on; ylabel('Max UEs in 1 Beam');
xlabel('Time (minutes from start)');
title('Beam Congestion Level');

exportgraphics(figA, fullfile(out_dir, 'peak_sat_resource_load.png'), 'Resolution', FIG.dpi);

%% Waterfall: Spatial Congestion vs Time
figure('Color', 'w', 'Name', 'Beam Load Waterfall', 'Position', FIG.wide);
imagesc(time_mins, 1:num_beams, beam_ue_count_matrix);
colormap(parula);
cb = colorbar;
ylabel(cb, 'UEs per Beam');
xlabel('Time (mins)'); ylabel('Beam ID');
title(sprintf('Sat %d: Spatial Congestion Waterfall', peakSat));

exportgraphics(gcf, fullfile(out_dir, 'peak_sat_congestion_waterfall.png'), 'Resolution', FIG.dpi);

%% Spatial Load Heatmap (Steering Angle Perspective)
counts_at_peak = beam_ue_count_matrix(:, t_peak);
active_beam_indices = find(counts_at_peak > 0);
active_counts = counts_at_peak(active_beam_indices);

if isempty(active_counts)
    cap_val = 1;
else
    cap_val = prctile(active_counts, 90);
    if cap_val == 0; cap_val = max(active_counts); end
end
fprintf('Capping color scale at %d UEs (90th Percentile).\n', round(cap_val));

b_eta = asind(sqrt(b_u.^2 + b_v.^2));
b_az  = atan2d(b_v, b_u);
b_x_deg = b_eta .* cosd(b_az);
b_y_deg = b_eta .* sind(b_az);
r_beam_deg = Beamwidth_deg / 2;

theta = linspace(0, 2*pi, 40);
X_all = zeros(length(theta), num_beams);
Y_all = zeros(length(theta), num_beams);

for b = 1:num_beams
    eta = b_eta(b);
    az  = b_az(b);
    if eta < 1e-3
        map_distortion_factor = 1;
    else
        map_distortion_factor = (eta * pi/180) / sind(eta);
    end
    r_major = r_beam_deg / cosd(eta);
    r_minor = r_beam_deg * map_distortion_factor;
    x_local = r_major * cos(theta);
    y_local = r_minor * sin(theta);
    x_rot = x_local * cosd(az) - y_local * sind(az);
    y_rot = x_local * sind(az) + y_local * cosd(az);
    X_all(:, b) = b_x_deg(b) + x_rot;
    Y_all(:, b) = b_y_deg(b) + y_rot;
end

X_active = X_all(:, active_beam_indices);
Y_active = Y_all(:, active_beam_indices);

% PLOT A: FULL FOV SPATIAL LOAD
fig_spatial = figure('Color', 'w', 'Visible', 'off', 'Position', FIG.square);
hold on; axis equal; box on;

plot(X_all, Y_all, 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5);
patch(X_active, Y_active, active_counts', 'EdgeColor', 'k', 'LineWidth', 0.5, 'FaceAlpha', 0.85);
plot(eta_max*cos(theta), eta_max*sin(theta), 'k--', 'LineWidth', 2);

colormap(turbo);
try clim([0 cap_val]); catch; caxis([0 cap_val]); end
cb = colorbar;
ylabel(cb, sprintf('Connected UEs (Capped at 90th pct: %d)', round(cap_val)));
xlabel('X Steering Angle (Degrees off Nadir)');
ylabel('Y Steering Angle (Degrees off Nadir)');
title(sprintf('Sat %d: FoV Beam Load Heatmap\nTime: %s', peakSat, datestr(time_vec(t_peak))));
xlim([-eta_max-2, eta_max+2]); ylim([-eta_max-2, eta_max+2]);

exportgraphics(fig_spatial, fullfile(out_dir, 'peak_sat_spatial_full_deg.png'), 'Resolution', FIG.dpi);

% PLOT B: ZOOMED "FREQUENCY REUSE" VIEW
[max_ue_val, max_idx] = max(active_counts);
b_idx_max = active_beam_indices(max_idx);
center_x_deg = b_x_deg(b_idx_max);
center_y_deg = b_y_deg(b_idx_max);

fprintf('Zooming in on the busiest beam: ID %d with %d UEs.\n', b_idx_max, max_ue_val);

fig_zoom = figure('Color', 'w', 'Visible', 'off', 'Position', FIG.square);
hold on; axis equal; box on; grid on;

plot(X_all, Y_all, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
patch(X_active, Y_active, active_counts', 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.85);

for i = 1:length(active_beam_indices)
    b_idx = active_beam_indices(i);
    text(b_x_deg(b_idx), b_y_deg(b_idx), sprintf('ID: %d\nUEs: %d', b_idx, counts_at_peak(b_idx)), ...
        'HorizontalAlignment', 'center', 'FontSize', FIG.annot_size, 'FontWeight', 'bold', 'Color', 'w', 'Clipping', 'on');
end

eta_max_beam = b_eta(b_idx_max);
az_max_beam  = b_az(b_idx_max);
if eta_max_beam < 1e-3
    map_dist_max = 1;
else
    map_dist_max = (eta_max_beam * pi/180) / sind(eta_max_beam);
end
r_major_max = r_beam_deg / cosd(eta_max_beam);
r_minor_max = r_beam_deg * map_dist_max;
x_loc_max = r_major_max * cos(theta);
y_loc_max = r_minor_max * sin(theta);
x_rot_max = x_loc_max * cosd(az_max_beam) - y_loc_max * sind(az_max_beam);
y_rot_max = x_loc_max * sind(az_max_beam) + y_loc_max * cosd(az_max_beam);
plot(center_x_deg + x_rot_max, center_y_deg + y_rot_max, 'r', 'LineWidth', 3);

colormap(turbo);
try clim([0 cap_val]); catch; caxis([0 cap_val]); end
cb_zoom = colorbar; ylabel(cb_zoom, 'Connected UEs');

zoom_radius = 2 * (r_beam_deg * 2);
xlim([center_x_deg - zoom_radius, center_x_deg + zoom_radius]);
ylim([center_y_deg - zoom_radius, center_y_deg + zoom_radius]);
xlabel('X Steering Angle (deg)'); ylabel('Y Steering Angle (deg)');
title(sprintf('Zoomed Active Region: Sat %d\nCentered on Peak Beam %d (Max Load)', peakSat, b_idx_max));

exportgraphics(fig_zoom, fullfile(out_dir, 'peak_sat_spatial_zoomed_deg.png'), 'Resolution', FIG.dpi);

%% Beam Load Percentile Distribution
counts_at_peak = beam_ue_count_matrix(:, t_peak);
active_counts = counts_at_peak(counts_at_peak > 0);
sorted_counts = sort(active_counts, 'ascend');

num_active = length(sorted_counts);
if num_active > 1
    percentiles = linspace(0, 100, num_active);
else
    percentiles = 100;
end
total_ues_connected = sum(sorted_counts);

fig_hist = figure('Color', 'w', 'Visible', 'off', 'Position', FIG.compact);
area(percentiles, sorted_counts, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', [0 0.3 0.6], 'LineWidth', 1.5);
hold on; grid on;

mean_load = mean(sorted_counts);
p90_load = prctile(sorted_counts, 90);
yline(mean_load, '--r', sprintf('Mean Load (%.1f UEs)', mean_load), 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
yline(p90_load, '--k', sprintf('90th Percentile (%.1f UEs)', p90_load), 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');

xlabel('Active Beam Percentile (%)');
ylabel('Number of Connected UEs');
title(sprintf('Sat %d: Beam Load Distribution\nTime: %s | Total Connected UEs: %d (Across %d Beams)', ...
    peakSat, datestr(time_vec(t_peak)), total_ues_connected, num_active));
xlim([0 100]);
ylim([0 max(sorted_counts) * 1.1]);

exportgraphics(fig_hist, fullfile(out_dir, 'peak_sat_load_percentile.png'), 'Resolution', FIG.dpi);
fprintf('Percentile Load Distribution plot saved.\n');


%% Beam Load CDF
counts_at_peak = beam_ue_count_matrix(:, t_peak);
active_counts = counts_at_peak(counts_at_peak > 0);
sorted_counts = sort(active_counts, 'ascend');

num_active = length(sorted_counts);
if num_active > 0
    percentiles = (1:num_active) / num_active * 100;
else
    percentiles = 0;
end
total_ues_connected = sum(sorted_counts);

fig_cdf = figure('Color', 'w', 'Visible', 'off', 'Position', FIG.tall);
plot(sorted_counts, percentiles, '-', 'LineWidth', 3, 'Color', [0 0.447 0.741]);
hold on; grid on;
area(sorted_counts, percentiles, 'FaceColor', [0 0.447 0.741], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

mean_load = mean(sorted_counts);
p90_load = prctile(sorted_counts, 90);
xline(mean_load, '--r', sprintf('Mean (%.1f UEs)', mean_load), ...
    'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', FIG.annot_size);
xline(p90_load, '--k', sprintf('90th Percentile (%.1f UEs)', p90_load), ...
    'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', FIG.annot_size);

xlabel('Number of Connected UEs per Beam');
ylabel('Cumulative Percentage of Beams (%)');
title(sprintf('Sat %d: Beam Load CDF\nTime: %s | Total Connected UEs: %d/%d (Across %d Beams)', ...
    peakSat, datestr(time_vec(t_peak)), total_ues_connected, Cfg.NumUEs, num_active));
if max(sorted_counts) > 0
    xlim([0 max(sorted_counts) * 1.05]);
else
    xlim([0 1]);
end
ylim([0 100]);

exportgraphics(fig_cdf, fullfile(out_dir, 'peak_sat_load_cdf.png'), 'Resolution', FIG.dpi);
fprintf('CDF Load Distribution plot saved.\n');
fprintf('All figures saved to: %s\n', out_dir);
