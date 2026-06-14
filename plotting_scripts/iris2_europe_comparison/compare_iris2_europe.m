% COMPARE_IRIS2_EUROPE  Coverage comparison of 4 IRIS2 constellation candidates over Europe.
%   Runs Constellation_simulator for each of the four IRIS2 candidate constellations
%   (Global LEO Star, Global LEO Delta, Regional LEO Delta, Global MEO) using a
%   shared population-based UE distribution over Europe (20 000 people per UE).
%   Calls plot_simulation per constellation, then saves a coverage summary figure.
%
% Headless run:
%   xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nodesktop -batch "compare_iris2_europe"

clear all; close all; clc;
delete(gcp('nocreate'));

% Add the project to the MATLAB path (robust to the script's folder depth).
repo_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
addpath(fullfile(repo_root, 'functions'));
path_setup();

%% ===== SHARED PARAMETERS =====
Lat_range_deg = [34, 72];
Lon_range_deg = [-12, 45];

StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
StopTime  = datetime('2-Jun-2025 14:59:59', 'TimeZone', 'UTC');

% Ku-band link budget (shared across all constellations)
f_DL           = 12e9;
B_DL           = 250e6;
NF_DL          = 5;
G_rx           = 33;        % dBi, UE receive antenna gain
Target_PFD_MHz = -120;      % dBW/m²/MHz
Target_Beamsize_Radius_km = 30;

FRF = 3;
RU  = 1;

people_per_ue = 20000;

use_parallel = false;
calc_link    = true;

%% ===== CONSTELLATION DEFINITIONS =====
% Option Num_planes Sats_per_plane Inclination Phasing Total_sats Height_km  MinElev
constellations(1).label          = 'Global LEO Star';
constellations(1).Num_planes     = 12;
constellations(1).Sats_per_plane = 22;
constellations(1).Inclination    = 90;
constellations(1).Phasing        = 6;
constellations(1).Total_sats     = 264;
constellations(1).height_km      = 1200;
% constellations(1).min_elev       = 38;
constellations(1).min_elev       = 20;
constellations(1).WalkerStar     = true;
constellations(1).dir_name       = 'global_leo_star';

constellations(2).label          = 'Global LEO Delta';
constellations(2).Num_planes     = 24;
constellations(2).Sats_per_plane = 11;
constellations(2).Inclination    = 75;
constellations(2).Phasing        = 10;
constellations(2).Total_sats     = 264;
constellations(2).height_km      = 1200;
% constellations(2).min_elev       = 30;
constellations(2).min_elev       = 20;
constellations(2).WalkerStar     = false;
constellations(2).dir_name       = 'global_leo_delta';

constellations(3).label          = 'Regional LEO Delta';
constellations(3).Num_planes     = 22;
constellations(3).Sats_per_plane = 12;
constellations(3).Inclination    = 58.5;
constellations(3).Phasing        = 9;
constellations(3).Total_sats     = 264;
constellations(3).height_km      = 1200;
% constellations(3).min_elev       = 42;
constellations(3).min_elev       = 20;
constellations(3).WalkerStar     = false;
constellations(3).dir_name       = 'regional_leo_delta';

constellations(4).label          = 'Global MEO';
constellations(4).Num_planes     = 6;
constellations(4).Sats_per_plane = 3;
constellations(4).Inclination    = 55;
constellations(4).Phasing        = 2;
constellations(4).Total_sats     = 18;
constellations(4).height_km      = 8000;
% constellations(4).min_elev       = 28;
constellations(4).min_elev       = 20;
constellations(4).WalkerStar     = false;
constellations(4).dir_name       = 'global_meo';

nConstellations = numel(constellations);

%% ===== OUTPUT DIRECTORY =====
script_dir = fileparts(mfilename('fullpath'));
% date_str   = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
run_dir    = fullfile(script_dir, '..', 'figures', 'iris2_europe');
for ci = 1:nConstellations
    mkdir(fullfile(run_dir, constellations(ci).dir_name));
end
fprintf('Saving all outputs to: %s\n', run_dir);

%% ===== GENERATE POPULATION-BASED UEs (shared across all constellations) =====
fprintf('\nGenerating population-based UEs over Europe (%d people/UE)...\n', people_per_ue);
[UE_lats, UE_lons, Total_Pop] = generate_population_ues(Lat_range_deg, Lon_range_deg, people_per_ue);
fprintf('Generated %d UEs representing %.2e people.\n', numel(UE_lats), Total_Pop);

%% ===== RUN SIMULATIONS =====
metrics_all = cell(nConstellations, 1);

for ci = 1:nConstellations
    con = constellations(ci);
    fprintf('\n========== [%d/%d] %s ==========\n', ci, nConstellations, con.label);

    Cfg = build_cfg(con, UE_lats, UE_lons, StartTime, StopTime, ...
                    f_DL, B_DL, NF_DL, G_rx, Target_PFD_MHz, FRF, RU, Target_Beamsize_Radius_km);
    Cfg.Total_Pop = Total_Pop;
    Cfg.Save_dir  = fullfile(run_dir, con.dir_name);

    metrics = Constellation_simulator(Cfg, use_parallel, calc_link);
    metrics_all{ci} = metrics;
    fprintf('%s -> worst coverage: %.4f%%\n', con.label, metrics.worst_coverage_percent);

    plot_simulation(metrics, use_parallel);
end

%% ===== SUMMARY FIGURE =====
labels        = {constellations.label};
worst_cov     = cellfun(@(m) m.worst_coverage_percent,     metrics_all);
mean_cov      = cellfun(@(m) mean(m.prob_coverage, 'all'), metrics_all);
total_sats_v  = [constellations.Total_sats];

fig_sum = figure('Color', 'w', 'Position', [100 100 700 420]);

subplot(1,2,1);
b = bar([worst_cov; mean_cov]', 'grouped');
set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 20, 'FontSize', 10);
ylabel('Coverage (%)');
title('Coverage over Europe');
legend({'Worst-case', 'Mean'}, 'Location', 'southeast');
grid on; box on;
ylim([0, 100]);

subplot(1,2,2);
bar(worst_cov ./ total_sats_v);
set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 20, 'FontSize', 10);
ylabel('Worst-case coverage / satellite (%)');
title('Coverage efficiency');
grid on; box on;

sgtitle('IRIS2 Constellation Candidates — Europe', 'FontSize', 13, 'FontWeight', 'bold');

fig_out = fullfile(run_dir, 'IRIS2_Europe_Summary.png');
exportgraphics(fig_sum, fig_out, 'Resolution', 300);
fprintf('\nSummary figure saved: %s\n', fig_out);

%% ===== TEXT SUMMARY =====
fprintf('\n%-22s  %6s  %8s  %8s  %6s\n', 'Constellation', 'Sats', 'Worst(%)', 'Mean(%)', 'Cov/sat');
fprintf('%s\n', repmat('-', 1, 62));
for ci = 1:nConstellations
    con = constellations(ci);
    m   = metrics_all{ci};
    fprintf('%-22s  %6d  %8.3f  %8.3f  %6.3f\n', ...
        con.label, con.Total_sats, m.worst_coverage_percent, ...
        mean(m.prob_coverage, 'all'), m.worst_coverage_percent / con.Total_sats);
end
fprintf('\nDone. All outputs saved to: %s\n', run_dir);

%% =================================================================
function Cfg = build_cfg(con, UE_lats, UE_lons, StartTime, StopTime, ...
                          f_DL, B_DL, NF_DL, G_rx, Target_PFD_MHz, FRF, RU, Target_Beamsize_Radius_km)
    h_m = con.height_km * 1e3;

    Cfg.Orbit_height               = h_m;
    Cfg.Min_elevation_UE           = con.min_elev;
    Cfg.SampleTime                 = 660;
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

    G_tx = gain_from_nadir_beam_size(con.height_km, Target_Beamsize_Radius_km);

    Cfg.DL.Direction    = "DL";
    Cfg.DL.f            = f_DL;
    Cfg.DL.B            = B_DL;
    Cfg.DL.NF           = NF_DL;
    Cfg.DL.G_rx         = G_rx;
    Cfg.DL.Tx_type      = "array";
    Cfg.DL.Rx_type      = "array";
    Cfg.DL.G_tx         = G_tx;
    Cfg.DL.BeamGrid     = calculate_hexagonal_beams(G_tx, f_DL, h_m, con.min_elev, Target_PFD_MHz, FRF);
    Cfg.DL.Max_EIRP_dBm = max(Cfg.DL.BeamGrid.BeamCenter_EIRP_dBmHz(:)) + 60;

    Cfg.WalkerStar      = con.WalkerStar;
    Cfg.Num_planes      = con.Num_planes;
    Cfg.Sats_per_plane  = con.Sats_per_plane;
    Cfg.Total_sats      = con.Total_sats;
    Cfg.Inclination     = con.Inclination;
    Cfg.Phasing         = con.Phasing;
end
