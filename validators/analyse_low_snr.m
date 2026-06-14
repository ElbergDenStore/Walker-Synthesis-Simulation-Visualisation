% ANALYSE_LOW_SNR  Diagnose where and why downlink SNR is low across the coverage area.
%
% Decomposes per-(UE, time) SNR into its additive link-budget terms to explain
% the gap between min and mean SNR. The identity that must hold exactly:
%
%   SNR = EIRP_center + beam_pointing_loss + G_rx − FSPL − Absorption − Rx_steering − Noise
%
% Because EIRP_center is scaled to a fixed PFD target, the (EIRP_center − FSPL)
% term is ~constant across the footprint. The two position-varying penalties are:
%   * beam_pointing_loss : AF+EF roll-off from serving beam centre (beam placement)
%   * Rx_steering_loss   : uncompensated UE-array scan loss at low elevation
%
% Outputs:
%   simulation_output/LowSNR_<timestamp>/LowSNR_drivers.png
%   simulation_output/LowSNR_<timestamp>/LowSNR_map.png
%
% Run:
%   Press play (F5) in VS Code, or:
%   xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nodesktop -batch "run('validators/analyse_low_snr.m')"

clear; close all; clc;

repo_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
addpath(fullfile(repo_root, 'functions'));
path_setup();

%% ===== CONFIGURE HERE =====
height_km        = 1000;
min_elevation_UE = 20;

Lat_range_deg = [54+(35/60), 83+(40/60)];
Lon_range_deg = [-(73+(10/60)), 33+(30/60)];

StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
StopTime  = datetime('1-Jun-2025 14:59:59', 'TimeZone', 'UTC');

f_DL           = 12e9;
B_DL           = 250e6;
NF_DL          = 5;
G_rx           = 33;
Target_PFD_MHz = -125;
FRF            = 3;
RU             = 1;

NumUEs = 100;

% Walker Star
Num_planes     = 5;
Sats_per_plane = 13;
Inclination    = 90;

low_pct = 5;   % "low SNR" = bottom N% of all (UE,time) samples

%% ===== BUILD AND RUN SIMULATION =====
[UE_lats, UE_lons] = generate_equal_area_ues(Lat_range_deg, Lon_range_deg, NumUEs);

Cfg.Orbit_height               = height_km * 1e3;
Cfg.Min_elevation_UE           = min_elevation_UE;
Cfg.SampleTime                 = 60;
Cfg.FRF                        = FRF;
Cfg.RU                         = RU;
Cfg.Use_P618                   = false;
Cfg.Modified_shannon           = true;
Cfg.Simple_Atmospheric_Loss_dB = 1;
Cfg.Share_bandwidth            = true;
Cfg.Target_PFD_MHz             = Target_PFD_MHz;
Cfg.StartTime                  = StartTime;
Cfg.StopTime                   = StopTime;
Cfg.Flat_UE_array.Lats         = UE_lats;
Cfg.Flat_UE_array.Lons         = UE_lons;

Cfg.DL.Direction    = "DL";
Cfg.DL.f            = f_DL;
Cfg.DL.B            = B_DL;
Cfg.DL.NF           = NF_DL;
Cfg.DL.G_rx         = G_rx;
Cfg.DL.Tx_type      = "array";
Cfg.DL.Rx_type      = "array";
Cfg.DL.G_tx         = get_adjusted_tx_gain(Cfg.Orbit_height, min_elevation_UE, f_DL);
Cfg.DL.BeamGrid     = calculate_hexagonal_beams(Cfg.DL.G_tx, f_DL, Cfg.Orbit_height, min_elevation_UE, Target_PFD_MHz, FRF);
Cfg.DL.Max_EIRP_dBm = max(Cfg.DL.BeamGrid.BeamCenter_EIRP_dBmHz(:));

Cfg.WalkerStar      = true;
Cfg.Num_planes      = Num_planes;
Cfg.Sats_per_plane  = Sats_per_plane;
Cfg.Total_sats      = Num_planes * Sats_per_plane;
Cfg.Inclination     = Inclination;
Cfg.Phasing         = Num_planes / 2;

fprintf('Running simulation (%d UEs, %d sats, %d km)...\n', NumUEs, Cfg.Total_sats, height_km);
metrics = Constellation_simulator(Cfg, false, true);

%% ===== FLATTEN PER-(UE, TIME) MATRICES =====
Re = 6378.137e3;
BG = Cfg.DL.BeamGrid;
UEs = metrics.UEs;

SD  = [UEs.SimData];
el  = vertcat(SD.Elevation_deg);
az  = vertcat(SD.Azimuth_deg);

DL    = [UEs.DL];
SNR   = vertcat(DL.SNR);
SINR  = vertcat(DL.SINR);
carr  = vertcat(DL.Carrier_density_dBmHz);
noise = vertcat(DL.Noise_density_dBmHz);
intf  = vertcat(DL.Interference_density_dBmHz);
fspl  = vertcat(DL.FSPL);
tloss = vertcat(DL.Total_loss);
sbi   = vertcat(DL.serving_beam_idx);
sbs   = vertcat(DL.serving_beam_signal_lin);
Abs   = [DL.Absorption];
absat = vertcat(Abs.At);

[U, T] = size(SNR);
lat = [UEs.Lat]';
lon = [UEs.Lon]';

%% ===== BEAM-OFFSET GEOMETRY =====
% Same direction-cosine transform used inside beam_gain_and_interference.m
eta = asind((Re / (Re + Cfg.Orbit_height)) .* cosd(el));
u   = sind(eta) .* cosd(az + 180);
v   = sind(eta) .* sind(az + 180);

has_beam = isfinite(sbi) & sbi > 0;
uc = nan(size(sbi));  vc = uc;  eirp_c = uc;
uc(has_beam)     = BG.u_center(sbi(has_beam));
vc(has_beam)     = BG.v_center(sbi(has_beam));
eirp_c(has_beam) = BG.BeamCenter_EIRP_dBmHz(sbi(has_beam));

rho_off     = hypot(u - uc, v - vc);
off_ang_deg = asind(min(rho_off, 1));

% Reference radii (for plot lines)
r_beam_uv = 1.391 / (BG.Nu * (pi/2));
r_pack_uv = sqrt(3) * r_beam_uv / 2;

%% ===== ADDITIVE SNR DECOMPOSITION =====
beam_loss    = sbs - eirp_c;        % AF+EF pointing loss (dB, <= 0)
rx_steer     = tloss - fspl - absat; % UE scan loss (dB, >= 0)
eirp_m_fspl  = eirp_c - fspl;       % PFD-compensated term (~constant)

% Verify reconstruction (residual should be ~0).
snr_model = eirp_c + beam_loss + G_rx - fspl - absat - rx_steer - noise;
recon_err = SNR - snr_model;

%% ===== LOW-SNR TAIL =====
valid      = isfinite(SNR);
snr_valid  = SNR(valid);
thr        = prctile(snr_valid, low_pct);
low        = valid & (SNR <= thr);
hi         = valid & (SNR >  thr);

%% ===== PER-TERM GROUP MEANS =====
term_names = {'Elevation_deg','BeamOffset_deg','BeamPointingLoss', ...
              'RxSteeringLoss','Absorption','EIRPc_minus_FSPL','Noise_dBmHz','Interf_dBmHz','SINR'};
term_mats  = {el, off_ang_deg, beam_loss, rx_steer, absat, eirp_m_fspl, noise, intf, SINR};
nTerms     = numel(term_names);
mean_low   = zeros(nTerms, 1);
mean_hi    = zeros(nTerms, 1);
for k = 1:nTerms
    M = term_mats{k};
    mean_low(k) = mean(M(low),  'omitnan');
    mean_hi(k)  = mean(M(hi),   'omitnan');
end
delta = mean_low - mean_hi;

%% ===== CORRELATIONS =====
corr_el  = local_corr(SNR(valid), el(valid));
corr_off = local_corr(SNR(valid), off_ang_deg(valid));
corr_bl  = local_corr(SNR(valid), beam_loss(valid));
corr_rxs = local_corr(SNR(valid), rx_steer(valid));

%% ===== PER-UE SUMMARY TABLE =====
snr_ue_min  = min(SNR, [], 2, 'omitnan');
snr_ue_mean = mean(SNR, 2, 'omitnan');
el_ue_mean  = mean(el, 2, 'omitnan');
off_ue_mean = mean(off_ang_deg, 2, 'omitnan');
UE_table = table((1:U)', lat, lon, snr_ue_min, snr_ue_mean, el_ue_mean, off_ue_mean, ...
    'VariableNames', {'UE','Lat','Lon','SNR_min','SNR_mean','Elev_mean','BeamOff_mean'});
UE_table = sortrows(UE_table, 'SNR_min', 'ascend');

%% ===== CONSOLE REPORT =====
fprintf('\n================ LOW-SNR DIAGNOSIS ================\n');
fprintf('Samples : %d UEs x %d steps = %d  (valid links: %d)\n', U, T, U*T, nnz(valid));
fprintf('SNR [dB]: min %.2f | 10%% %.2f | mean %.2f | median %.2f | max %.2f\n', ...
    min(snr_valid), prctile(snr_valid,10), mean(snr_valid), median(snr_valid), max(snr_valid));
fprintf('Reconstruction residual (should be ~0): max |err| = %.4f dB\n', max(abs(recon_err(valid))));
fprintf('Low-SNR tail = bottom %g%% (SNR <= %.2f dB), %d samples\n', low_pct, thr, nnz(low));

fprintf('\nPearson corr of SNR with:\n');
fprintf('   elevation          : %+.3f\n', corr_el);
fprintf('   beam offset angle  : %+.3f\n', corr_off);
fprintf('   beam pointing loss : %+.3f  (near +1 => beam PLACEMENT drives low SNR)\n', corr_bl);
fprintf('   Rx steering loss   : %+.3f  (near -1 => ELEVATION drives low SNR)\n', corr_rxs);

fprintf('\n%-22s  %9s  %9s  %9s\n', 'Term', 'low-SNR', 'rest', 'delta');
fprintf('%s\n', repmat('-', 1, 55));
for k = 1:nTerms
    fprintf('%-22s  %9.2f  %9.2f  %+9.2f\n', term_names{k}, mean_low(k), mean_hi(k), delta(k));
end
fprintf('\nWorst 10 UEs by min SNR:\n');
disp(UE_table(1:min(10,U), :));

%% ===== OUTPUT DIRECTORY =====
date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir  = fullfile(repo_root, 'simulation_output', sprintf('LowSNR_%s', date_str));
mkdir(out_dir);
fprintf('Saving figures to: %s\n', out_dir);

%% ===== FIGURE 1: SNR drivers =====
% Pre-extract valid samples as flat vectors so subsample_mask indexing works cleanly.
el_v   = el(valid);
snr_v  = SNR(valid);
off_v  = off_ang_deg(valid);
sel    = subsample_mask(numel(el_v), 20000);

f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 800]);
tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) SNR vs elevation
nexttile;
scatter(el_v(sel), snr_v(sel), 6, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.15);
grid on; box on; xlabel('Elevation (deg)'); ylabel('SNR (dB)');
hold on;
[med_x, med_y] = binned_median(el_v, snr_v, 20);
plot(med_x, med_y, 'r-', 'LineWidth', 2);
title(sprintf('SNR vs elevation (\\rho = %+.2f)', corr_el));

% (b) SNR vs beam offset
nexttile;
scatter(off_v(sel), snr_v(sel), 6, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.15);
grid on; box on; xlabel('Offset from serving beam centre (deg)'); ylabel('SNR (dB)');
hold on;
[med_x, med_y] = binned_median(off_v, snr_v, 20);
plot(med_x, med_y, 'r-', 'LineWidth', 2);
xline(asind(min(r_beam_uv,1)),  'r--', '3 dB radius',     'LabelVerticalAlignment','bottom');
xline(asind(min(r_pack_uv,1)),  'm--', 'hex worst-case',  'LabelVerticalAlignment','bottom');
title(sprintf('SNR vs beam placement (\\rho = %+.2f)', corr_off));

% (c) Beam-offset distribution
nexttile;
histogram(off_v, 50, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
hold on;
xline(asind(min(r_beam_uv,1)), 'r--', '3 dB radius');
xline(asind(min(r_pack_uv,1)), 'm--', 'hex worst-case');
grid on; box on;
xlabel('Offset from serving beam centre (deg)'); ylabel('Sample count');
title('Beam-offset distribution (gaps \Rightarrow coverage holes)');

% (d) Component delta bar chart: what differs in the low-SNR tail
nexttile;
comp_names = {'BeamPointingLoss','RxSteeringLoss','Absorption','EIRPc_minus_FSPL','Noise_dBmHz'};
comp_idx   = cellfun(@(s) find(strcmp(term_names, s)), comp_names);
barh(delta(comp_idx));
set(gca, 'YTick', 1:numel(comp_names), 'YTickLabel', strrep(comp_names, '_', '\_'));
grid on; box on;
xlabel('low-SNR mean \minus rest mean (dB)');
title('What differs in the low-SNR tail');

sgtitle(sprintf('Low-SNR drivers  |  %d sats, %.0f km, min elev %.0f\\circ, FRF %d', ...
    Cfg.Total_sats, height_km, Cfg.Min_elevation_UE, Cfg.FRF), 'FontWeight', 'bold');
exportgraphics(f1, fullfile(out_dir, 'LowSNR_drivers.png'), 'Resolution', 200);
close(f1);

%% ===== FIGURE 2: Geographic map of min SNR per UE =====
f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 700]);
hold on;
try
    land = shaperead('landareas.shp', 'UseGeoCoords', true);
    for i = 1:numel(land)
        plot(land(i).Lon, land(i).Lat, 'Color', [0.6 0.6 0.6]);
    end
catch
end
scatter(lon, lat, 60, snr_ue_min, 'filled', 'MarkerEdgeColor', 'k');
nWorst = min(5, U);
plot(UE_table.Lon(1:nWorst), UE_table.Lat(1:nWorst), 'rp', ...
    'MarkerSize', 16, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
cb = colorbar; ylabel(cb, 'Worst-case (min) SNR per UE (dB)');
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
xlim([min(lon)-2, max(lon)+2]); ylim([min(lat)-2, max(lat)+2]);
grid on; box on;
title('Where low SNR happens (red stars = 5 worst UEs)');
exportgraphics(f2, fullfile(out_dir, 'LowSNR_map.png'), 'Resolution', 200);
close(f2);

fprintf('Done.\n');

%% ===== LOCAL HELPERS (inlined — scripts cannot have nested functions) =====

function r = local_corr(x, y)
    m = isfinite(x) & isfinite(y);
    if nnz(m) < 2; r = NaN; return; end
    c = corrcoef(x(m), y(m));
    r = c(1, 2);
end

function mask = subsample_mask(n_total, n_max)
    mask = false(n_total, 1);
    if n_total <= n_max
        mask(:) = true;
    else
        mask(randperm(n_total, n_max)) = true;
    end
end

function [ctrs, med] = binned_median(x, y, nbins)
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    edges = linspace(min(x), max(x), nbins + 1);
    ctrs  = movmean(edges, 2, 'Endpoints', 'discard');
    med   = nan(size(ctrs));
    for k = 1:numel(ctrs)
        in = x >= edges(k) & x < edges(k+1);
        if any(in), med(k) = median(y(in)); end
    end
end
