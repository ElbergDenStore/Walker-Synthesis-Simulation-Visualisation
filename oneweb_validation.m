
% oneweb_validation
% Standalone OneWeb-style validation scenario without get_cfg.m.
% Uses 55 satellites/plane for simplicity (12 planes => 660 satellites).

clearvars; close all; clc;
% path_setup()

%% 1) OneWeb-like constellation and simulation settings
Cfg = struct();
Cfg.WalkerStar = true;
Cfg.Orbit_height = 1200e3;      % meters
Cfg.Inclination = 87.9;         % degrees
Cfg.Num_planes = 12;
Cfg.Sats_per_plane = 55;        % user choice: 55 per orbit
Cfg.Total_sats = Cfg.Num_planes * Cfg.Sats_per_plane;
Cfg.Phasing = floor(Cfg.Num_planes / 2);

Cfg.Min_elevation_UE = 20;
Cfg.SampleTime = 60;            % seconds (adjust if you need finer time resolution)
Cfg.StartTime = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');
Cfg.StopTime = datetime('1-Jun-2025 13:59:59', 'TimeZone', 'UTC'); %% TODO make the simulation time longer
Cfg.Modified_shannon = true;
Cfg.Share_bandwidth = true;
Cfg.Use_P618 = true;
Cfg.Ignore_Interference = true;

%% 2) Region and population-based UE generation
% Denmark-only region
latlim = [54.5, 58.0];
lonlim = [8.0, 15.5];
people_per_ue = 30000;

try
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, total_pop] = generate_population_based_UEs(latlim, lonlim, people_per_ue);
    Cfg.NumUEs = numel(Cfg.Flat_UE_array.Lats);
catch ME
    warning('%s', sprintf('generate_population_based_UEs failed: %s. Falling back to single UE in Aalborg for smoke test.', ME.message));
    % Aalborg approximate coordinates
    Cfg.Flat_UE_array.Lats = 57.0488;
    Cfg.Flat_UE_array.Lons = 9.9217;
    total_pop = NaN;
    Cfg.NumUEs = 1;
end

fprintf('UE population model: 1 UE per %d people\n', people_per_ue);
fprintf('Total population in bounds: %.0f\n', total_pop);
fprintf('Generated UEs: %d\n', Cfg.NumUEs);

%% 3) Link setup (Ku downlink center frequency 11.7 GHz)
Cfg.DL = struct();
Cfg.DL.Direction = "DL";
Cfg.DL.f = 11.7e9;
Cfg.DL.B = 250e6;                % 2GHz is possible, but does not seem to be used %https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10551683&tag=1
Cfg.DL.Rx_type = "array";

% Downlink EIRP density: -13.4 dBW / 4 kHz
% Convert to total EIRP over configured bandwidth B.
oneweb_eirp_density_dBW_per_4kHz = -13.4;
Cfg.DL.Max_EIRP_dBm = (oneweb_eirp_density_dBW_per_4kHz + 10*log10(Cfg.DL.B / 4e3)) + 30;
Cfg.DL.Max_EIRP_dBm_Hz = oneweb_eirp_density_dBW_per_4kHz + 30 - 10*log10(4e3);

% User terminal: Hughes HL1120W (dual 60x40 cm phased arrays), G/T up to 11.3 dB/K.
% Derive receive gain from physical aperture at 11.7 GHz, then back out NF from G/T. needed as Noise figure instead to get something realistic, instead of G/T which does not really mean anything
c = physconst('LightSpeed');
lambda = c / Cfg.DL.f;

eta_ap = 0.65;                       % assumed aperture efficiency
A_one_panel = 0.60 * 0.40;           % m^2
A_effective = 2 * A_one_panel;       % dual arrays
Cfg.DL.G_rx = 10*log10(eta_ap * 4*pi*A_effective / lambda^2);

GT_dBK = 11.3;                       % dB/K (Datasheet value)
Tsys_K = 10.^((Cfg.DL.G_rx - GT_dBK)/10);
Tant_K = 40;                         % clear-sky-ish antenna temp assumption
F_lin = 1 + max(Tsys_K - Tant_K, 1) / 290;
Cfg.DL.NF = 10*log10(F_lin);

% Keep a matching Tx gain field for beam geometry and steering loss internals.
Cfg.DL.BeamGrid = calculate_OneWeb_beams();

fprintf('OneWeb-style DL config:\n');
fprintf('  f = %.2f GHz, B = %.1f MHz\n', Cfg.DL.f/1e9, Cfg.DL.B/1e6);
fprintf('  Max EIRP = %.2f dBm (from %.1f dBW/4kHz)\n', Cfg.DL.Max_EIRP_dBm, oneweb_eirp_density_dBW_per_4kHz);
fprintf('  G_rx = %.2f dBi, NF = %.2f dB (from G/T = %.1f dB/K)\n', Cfg.DL.G_rx, Cfg.DL.NF, GT_dBK);

%% 4) Run simulation
calc_link = true;
use_parallel = true;
metrics = coverage_simulator_function(Cfg, use_parallel, calc_link);

%% 5) Basic output and optional plotting
fprintf('Simulation finished. Worst coverage: %.2f%%\n', metrics.worst_coverage_percent);
plot_simulation(metrics, use_parallel);



