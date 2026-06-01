% xvfb-run -a --server-args="-screen 0 1920x1080x24" matlab -nosplash -nodesktop -batch "run('plotting_scripts/plot_payload_power_model.m')"
close all; clearvars; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

%% Output directory
out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'plotting_scripts/figures/power');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% =========================================================================
% 5G NTN SATELLITE PAYLOAD POWER MODEL with element power cost
% =========================================================================

bandwidth = 50e6;
Linear_gain = linspace(10^(34*0.1), 10^(42*0.1), 500);
G_tx_dB = 10*log10(Linear_gain);
orbit_height = 1000e3;
f = 12e9;
Min_Elev_deg = 20;
FRF = 1;
Spectral_efficiency = 0.5; % bits/Hz

% --- DECOUPLED RF VARIABLES ---
Target_PFD_dBW_m2_MHz = -125;    % dBW/m²/MHz  (ITU regulatory limit for LEO)
RF_PA_efficiency = 0.10;

% Geometry — dimensioning on the worst-case UE: edge of coverage (min elevation)
Re_m       = 6371e3;
theta_dim  = asin(Re_m * cosd(Min_Elev_deg) / (Re_m + orbit_height));  % nadir angle at min elev
R_dim      = -Re_m*sind(Min_Elev_deg) + sqrt(Re_m^2*sind(Min_Elev_deg)^2 - Re_m^2 + (Re_m + orbit_height)^2);

% Spherical cap coverage area (for context)
rho_rad      = deg2rad(90 - Min_Elev_deg) - theta_dim;
coverage_km2 = 2 * pi * (Re_m/1e3)^2 * (1 - cos(rho_rad));

% PFD target in W/m² over the full bandwidth
PFD_target_lin = 10^(Target_PFD_dBW_m2_MHz / 10) * bandwidth / 1e6;

% Adaptive power control: for each beam position ρ, the satellite adjusts P_tx so
% the UE receives exactly PFD_target:
%   PFD_target = P_tx(ρ) × G_max × cos^1.5(θ(ρ)) / (4π R(ρ)²)
%   => P_tx(ρ) = PFD_target × 4π R(ρ)² / (G_max × cos^1.5(θ(ρ)))
%   => EIRP(ρ)  = PFD_target × 4π R(ρ)²  [varies with range, not gain]
%
% Average over the coverage footprint (uniform UE density on Earth surface):
rho_int   = linspace(0, rho_rad, 2000);
R_int     = sqrt(Re_m^2 + (Re_m + orbit_height)^2 - 2*Re_m*(Re_m + orbit_height)*cos(rho_int));
theta_int = asin(min(Re_m .* sin(rho_int) ./ R_int, 1));   % nadir angle for each ρ
w_int     = sin(rho_int);                                   % area weight ∝ sin(ρ) dρ

% Scalar factor: P_tx_avg = P_tx_avg_factor / G_max
P_tx_avg_factor = trapz(rho_int, (4*pi*PFD_target_lin .* R_int.^2 ./ cos(theta_int).^1.5) .* w_int) ...
                / trapz(rho_int, w_int);

% Average and peak EIRP (for reference; EIRP is independent of array gain)
EIRP_avg_linear = trapz(rho_int, (4*pi*PFD_target_lin .* R_int.^2) .* w_int) ...
                / trapz(rho_int, w_int);
EIRP_avg_dBW    = 10*log10(EIRP_avg_linear);
EIRP_peak_dBW   = 10*log10(PFD_target_lin * 4*pi * R_dim^2);   % edge-of-coverage (PA headroom)

fprintf('Edge-of-coverage: nadir angle = %.1f deg, slant range = %.0f km\n', ...
    rad2deg(theta_dim), R_dim/1e3);
fprintf('Average EIRP over footprint: %.1f dBW | Peak EIRP (edge): %.1f dBW\n', EIRP_avg_dBW, EIRP_peak_dBW);

% --- HARDWARE OVERHEAD VARIABLES ---
digital_transceiver_power = 13/16;  % ADRV9040 8T8R transceiver
receiver_LNA_power = 0.1;           % Watts per antenna element

fivegNTN_power = 12;                % 12 watts per gbps
Beam_utilization = 0.2;
GMAC_power = 1.25*1e-3;             % 1.25 mW per G-MAC
hybrid_analog_elements_per_digital = 1;

element_factor_dB = 6;
num_elements = ceil(Linear_gain ./ 10^(element_factor_dB*0.1));

Num_beams = zeros(1, length(Linear_gain));
fprintf('Computing beam counts for %d gain values...\n', length(Linear_gain));
for idx = 1:length(Linear_gain)
    BeamGrid = calculate_hexagonal_beams(G_tx_dB(idx), f, orbit_height, Min_Elev_deg, [], FRF);
    Num_beams(idx) = BeamGrid.num_beams;
end

% --- THROUGHPUT CALCULATIONS ---
Throughput = Spectral_efficiency * bandwidth * Num_beams * Beam_utilization;
Throughput_Gbps = Throughput ./ 1e9;
protocol_stack_power = fivegNTN_power * Throughput_Gbps;

digital_elements = ceil(num_elements ./ hybrid_analog_elements_per_digital);
digital_beamforming_power = GMAC_power .* digital_elements .* (2 * bandwidth / 1e9) .* Num_beams .* Beam_utilization;
processing_power = protocol_stack_power + digital_beamforming_power;

% --- SYSTEMATIC RF POWER CALCULATIONS (adaptive power control) ---
% P_tx varies per beam position; use area-weighted average over the footprint.
% P_RF_per_beam = P_tx_avg_factor / G_max  →  as N grows (G_max ∝ N), P_RF ∝ 1/N
P_RF_per_beam      = P_tx_avg_factor ./ Linear_gain;
P_RF_total_array   = P_RF_per_beam .* Num_beams * Beam_utilization;
RF_frontend_power  = P_RF_total_array ./ RF_PA_efficiency;

P_RF_per_element_mW = (P_RF_total_array ./ num_elements) .* 1000;

% --- HARDWARE OVERHEAD ---
total_LNA_power        = receiver_LNA_power .* num_elements;
total_transceiver_power = digital_transceiver_power .* digital_elements;
Hardware_Overhead_Power = total_LNA_power + total_transceiver_power;

% --- TOTAL POWER ---
Total_power = RF_frontend_power + Hardware_Overhead_Power + processing_power;

desired_db_ticks = [30, 33, 35, 38, 40, 42];

% =========================================================================
% FIGURE 1: PAYLOAD POWER BREAKDOWN (STACKED AREA)
% =========================================================================
f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 500, 400]);
ax1_f1 = axes(f1, 'Color', 'w');

Power_Matrix = [RF_frontend_power; Hardware_Overhead_Power; protocol_stack_power; digital_beamforming_power]';
colororder(ax1_f1, [[0, 0.447, 0.741]; [0.4660, 0.6740, 0.1880]; [0.850, 0.325, 0.098]; [0.929, 0.694, 0.125]]);
area(ax1_f1, num_elements, Power_Matrix, 'EdgeColor', 'none');
legend(ax1_f1, {'RF Transmit PA', 'Transceivers & LNAs', '5G NTN Protocol', 'Digital Beamforming'}, ...
    'Location', 'northwest', 'FontWeight', 'bold');
xlabel(ax1_f1, 'Antenna Array Elements', 'FontWeight', 'bold', 'FontSize', 14);
ylabel(ax1_f1, 'Total Power (W)',  'FontWeight', 'bold', 'FontSize', 14);
grid(ax1_f1, 'on'); ax1_f1.Box = 'off';

ax2_f1 = axes(f1, 'Position', ax1_f1.Position, 'Color', 'none', ...
              'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
              'XColor', 'k', 'YColor', 'k');
linkaxes([ax1_f1, ax2_f1], 'x');
xlim(ax1_f1, [min(num_elements), max(num_elements)]);
ax2_f1.XTick = ceil(10.^((desired_db_ticks - element_factor_dB) * 0.1));
ax2_f1.XTickLabel = string(desired_db_ticks) + " dBi";
ax2_f1.YTick = [];
title(ax2_f1, sprintf('Payload Processing Power Breakdown\nTarget PFD: %.0f dBW/m²/MHz | PA Eff: %.0f%% | Beam Util: %.0f%% | Analog/Dig Ratio: %d', ...
    Target_PFD_dBW_m2_MHz, RF_PA_efficiency*100, Beam_utilization*100, hybrid_analog_elements_per_digital), ...
    'FontSize', 12, 'FontWeight', 'bold');
ax1_f1.Position = [0.12, 0.12, 0.75, 0.68];
ax2_f1.Position = ax1_f1.Position;

exportgraphics(f1, fullfile(out_dir, 'payload_power_breakdown.png'), 'Resolution', 300);
close(f1); fprintf('Saved: payload_power_breakdown.png\n');

% =========================================================================
% FIGURE 2: POWER EFFICIENCY (W/Gbps) vs TX GAIN
% =========================================================================
f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [150, 150, 500, 400]);
ax1_f2 = axes(f2, 'Color', 'w');

Power_Efficiency = Total_power ./ Throughput_Gbps;
plot(ax1_f2, num_elements, Power_Efficiency, '-', 'LineWidth', 3, 'Color', [0, 0.447, 0.741]);
xlabel(ax1_f2, 'Array Elements',             'FontWeight', 'bold', 'FontSize', 14);
ylabel(ax1_f2, 'Power Efficiency (W/Gbps)',  'FontWeight', 'bold', 'FontSize', 14);
grid(ax1_f2, 'on'); ax1_f2.Box = 'off';

ax2_f2 = axes(f2, 'Position', ax1_f2.Position, 'Color', 'none', ...
              'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
              'XColor', 'k', 'YColor', 'k');
linkaxes([ax1_f2, ax2_f2], 'x');
xlim(ax1_f2, [min(num_elements), max(num_elements)]);
ax2_f2.XTick = ceil(10.^((desired_db_ticks - element_factor_dB) * 0.1));
ax2_f2.XTickLabel = string(desired_db_ticks) + " dBi";
ax2_f2.YTick = [];
title(ax2_f2, sprintf('Network Power Efficiency vs. Tx Gain\nBeam Util: %.0f%% | Analog/Dig Ratio: %d', ...
    Beam_utilization*100, hybrid_analog_elements_per_digital), 'FontSize', 12, 'FontWeight', 'bold');
ax1_f2.Position = [0.12, 0.12, 0.75, 0.68];
ax2_f2.Position = ax1_f2.Position;

exportgraphics(f2, fullfile(out_dir, 'payload_power_efficiency.png'), 'Resolution', 300);
close(f2); fprintf('Saved: payload_power_efficiency.png\n');

% =========================================================================
% FIGURE 3: POWER PER BEAM (W/Beam) vs TX GAIN
% =========================================================================
f3 = figure('Color', 'w', 'Visible', 'off', 'Position', [200, 200, 500, 400]);
ax1_f3 = axes(f3, 'Color', 'w');

Power_per_Beam = Total_power ./ Num_beams;
plot(ax1_f3, num_elements, Power_per_Beam, '-', 'LineWidth', 3, 'Color', [0.850, 0.325, 0.098]);
xlabel(ax1_f3, 'Array Elements',                    'FontWeight', 'bold', 'FontSize', 14);
ylabel(ax1_f3, 'Total DC Power per Beam (W/Beam)',  'FontWeight', 'bold', 'FontSize', 14);
grid(ax1_f3, 'on'); ax1_f3.Box = 'off';

ax2_f3 = axes(f3, 'Position', ax1_f3.Position, 'Color', 'none', ...
              'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
              'XColor', 'k', 'YColor', 'k');
linkaxes([ax1_f3, ax2_f3], 'x');
xlim(ax1_f3, [min(num_elements), max(num_elements)]);
ax2_f3.XTick = ceil(10.^((desired_db_ticks - element_factor_dB) * 0.1));
ax2_f3.XTickLabel = string(desired_db_ticks) + " dBi";
ax2_f3.YTick = [];
title(ax2_f3, sprintf('Energy Cost Per Active Beam vs. Tx Gain\nBeam Util: %.0f%% | Analog/Dig Ratio: %d', ...
    Beam_utilization*100, hybrid_analog_elements_per_digital), 'FontSize', 12, 'FontWeight', 'bold');
ax1_f3.Position = [0.12, 0.12, 0.75, 0.68];
ax2_f3.Position = ax1_f3.Position;

exportgraphics(f3, fullfile(out_dir, 'payload_power_per_beam.png'), 'Resolution', 300);
close(f3); fprintf('Saved: payload_power_per_beam.png\n');

% =========================================================================
% FIGURE 4: THE N^2 ADVANTAGE (Target EIRP & Element RF Power)
% =========================================================================
f4 = figure('Color', 'w', 'Visible', 'off', 'Position', [250, 250, 500, 400]);
ax1_f4 = axes(f4, 'Color', 'w');

yyaxis(ax1_f4, 'left');
plot(ax1_f4, num_elements, repmat(EIRP_avg_dBW, 1, length(num_elements)), '-', 'LineWidth', 3, 'Color', [0.4940, 0.1840, 0.5560]);
ylab = ylabel(ax1_f4, 'Average EIRP over Footprint (dBW)', 'FontWeight', 'bold', 'FontSize', 14);
ax1_f4.YColor = [0.4940, 0.1840, 0.5560];
ylim(ax1_f4, [EIRP_avg_dBW-2, EIRP_peak_dBW+2]);

yyaxis(ax1_f4, 'right');
plot(ax1_f4, num_elements, P_RF_per_element_mW, '--', 'LineWidth', 3, 'Color', [0.4660, 0.6740, 0.1880]);
ylabel(ax1_f4, 'Required RF Power per Element (mW)', 'FontWeight', 'bold', 'FontSize', 14);
ax1_f4.YColor = [0.4660, 0.6740, 0.1880];

xlabel(ax1_f4, 'Array Elements', 'FontWeight', 'bold', 'FontSize', 14);
grid(ax1_f4, 'on'); ax1_f4.Box = 'off';

ax2_f4 = axes(f4, 'Position', ax1_f4.Position, 'Color', 'none', ...
              'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
              'XColor', 'k', 'YColor', 'none');
linkaxes([ax1_f4, ax2_f4], 'x');
xlim(ax1_f4, [min(num_elements), max(num_elements)]);
ax2_f4.XTick = ceil(10.^((desired_db_ticks - element_factor_dB) * 0.1));
ax2_f4.XTickLabel = string(desired_db_ticks) + " dBi";
ax2_f4.YTick = [];
title(ax2_f4, sprintf('PFD Target: %.0f dBW/m²/MHz | \\epsilon_{min} = %.0f°, R = %.0f km\nBeam Util: %.0f%% | Analog/Dig Ratio: %d', ...
    Target_PFD_dBW_m2_MHz, Min_Elev_deg, R_dim/1e3, Beam_utilization*100, hybrid_analog_elements_per_digital), 'FontSize', 12, 'FontWeight', 'bold');
ax1_f4.Position = [0.12, 0.12, 0.75, 0.68];
ax2_f4.Position = ax1_f4.Position;

exportgraphics(f4, fullfile(out_dir, 'payload_n2_advantage.png'), 'Resolution', 300);
close(f4); fprintf('Saved: payload_n2_advantage.png\n');

% =========================================================================
% FIGURE 5: THROUGHPUT & BANDWIDTH (CONSTANT EIRP/MHz)
% =========================================================================
Var_Bandwidth      = repmat(bandwidth, 1, length(Linear_gain));
Var_Throughput_Gbps = (Spectral_efficiency .* Var_Bandwidth .* Num_beams .* Beam_utilization) ./ 1e9;

f5 = figure('Color', 'w', 'Visible', 'off', 'Position', [300, 300, 500, 400]);
ax1_f5 = axes(f5, 'Color', 'w');

yyaxis(ax1_f5, 'left');
plot(ax1_f5, num_elements, Var_Throughput_Gbps, '-', 'LineWidth', 3, 'Color', [0, 0.447, 0.741]);
ylabel(ax1_f5, 'Throughput (Gbps)', 'FontWeight', 'bold', 'FontSize', 14);
ax1_f5.YColor = [0, 0.447, 0.741];

yyaxis(ax1_f5, 'right');
plot(ax1_f5, num_elements, Var_Bandwidth ./ 1e6, '--', 'LineWidth', 3, 'Color', [0.850, 0.325, 0.098]);
ylabel(ax1_f5, 'Bandwidth (MHz)', 'FontWeight', 'bold', 'FontSize', 14);
ax1_f5.YColor = [0.850, 0.325, 0.098];
ylim(ax1_f5, [(bandwidth/1e6)-10, (bandwidth/1e6)+10]);

xlabel(ax1_f5, 'Array Elements', 'FontWeight', 'bold', 'FontSize', 14);
grid(ax1_f5, 'on'); ax1_f5.Box = 'off';

ax2_f5 = axes(f5, 'Position', ax1_f5.Position, 'Color', 'none', ...
              'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
              'XColor', 'k', 'YColor', 'none');
linkaxes([ax1_f5, ax2_f5], 'x');
xlim(ax1_f5, [min(num_elements), max(num_elements)]);
ax2_f5.XTick = ceil(10.^((desired_db_ticks - element_factor_dB) * 0.1));
ax2_f5.XTickLabel = string(desired_db_ticks) + " dBi";
ax2_f5.YTick = [];
title(ax2_f5, sprintf('Throughput vs Array Size at PFD Limit\nPFD: %.0f dBW/m²/MHz | Beam Util: %.0f%% | Ratio: %d', ...
    Target_PFD_dBW_m2_MHz, Beam_utilization*100, hybrid_analog_elements_per_digital), 'FontSize', 12, 'FontWeight', 'bold');
ax1_f5.Position = [0.12, 0.12, 0.75, 0.68];
ax2_f5.Position = ax1_f5.Position;

exportgraphics(f5, fullfile(out_dir, 'payload_throughput_bandwidth.png'), 'Resolution', 300);
close(f5); fprintf('Saved: payload_throughput_bandwidth.png\n');

fprintf('All figures saved to: %s\n', out_dir);
