% livescript_simulation_v2.m
% Combined simulation script merging:
% 1. Grid-based simulation loop & stats from livescript_simulation
% 2. Improved Link Budget Calculation (DL/UL structs) from UE_focus_simulation

clear all;
close all;
clc;

FORCE_RERUN = false; 

dataFile = 'BigSimulationData_v2.mat';
% dataFile = 'SmallSimulationData_v2.mat';

if exist(dataFile, 'file') && ~FORCE_RERUN 
    fprintf('Loading saved simulation data from %s...\n', dataFile);
    load(dataFile);
else
    fprintf('Running new simulation...\n');
    fprintf('Defining constellation...\n')
    sc = satelliteScenario;
    
    sc.StartTime = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
    % sc.StopTime  = datetime('1-Jun-2025 01:59:59', 'TimeZone', 'UTC');
    sc.StopTime  = datetime('1-Jun-2025 23:59:59', 'TimeZone', 'UTC');
    
    % Sample Time
    sc.SampleTime = 20;
    step = sc.SampleTime;
    
    simTimes = sc.StartTime:seconds(step):sc.StopTime; 
    simTimes.TimeZone = 'UTC';
    
    % Elevation angle
    minElevationUE = 20;
    
    % Constellation Parameters
    orbit_height = 700e3;
    r_earth = 6378.14e3;
    % Walker-Delta Constellation: 64 satellites, 8 planes, phasing = 3, i = 75°
    sat = walkerDelta(sc, orbit_height+r_earth, 75, 64, 8, 3, Name="S4D", OrbitPropagator="sgp4"); 
    
    %% Constellation Visualization
    % sensors = conicalSensor(sat, 'MaxViewAngle', 57); 
    % fieldOfView(sensors);
    % satelliteScenarioViewer(sc);
    
    %% UE/GS Grid
    lat_vec = linspace(50, 85, 16); 
    % lat_vec = linspace(50, 85, 1); 
    lon_vec = linspace(-60, 30, 3);
    [LonGrid, LatGrid] = meshgrid(lon_vec, lat_vec);
    nPts = numel(LonGrid);
    
    %% Coverage Simulation
    fprintf('Defining UEs and links to the satellite...\n')
    
    SatData = cell(1, nPts); 
    TimeDurationStats = cell(1, nPts); 
    
    dq = parallel.pool.DataQueue;
    fprintf('Simulation Progress:  0.00%%');
    afterEach(dq, @(~) updateCommandLine(nPts));
    
    tic
    parfor pt = 1:nPts
        lat = LatGrid(pt);
        lon = LonGrid(pt);
        ue_name_str = sprintf('UE%dLat%.0f', pt, lat);
        ue = groundStation(sc, lat, lon, 'Name', ue_name_str, 'MinElevationAngle', minElevationUE);
        
        ac = access(sat, ue);
        intvls = accessIntervals(ac);
        
        %%%%% Logic for TimeDurationStats
        if isempty(intvls)
            time_stats_struct = struct('Mean', 0, 'Max', 0, 'Min', 0, 'Total', 0, 'Num_Passes', 0);
        else
            durations = seconds(intvls.EndTime - intvls.StartTime);
            
            time_stats_struct = struct(...
                'Mean',   mean(durations), ...
                'Max',    max(durations), ...
                'Min',    min(durations), ...
                'Total',  sum(durations), ...
                'Num_Passes', numel(durations) ... 
            );
        end
        
        TimeDurationStats{pt} = time_stats_struct;
    
        %%%%% Logic for AER (azimuth, elevation and range)
        SatData{pt} = [];
        
        if ~isempty(intvls)
            % Estimate allocation size
            raw_steps = seconds(intvls.EndTime - intvls.StartTime) ./ sc.SampleTime;
            num_steps = round(raw_steps); 
            n_approx  = sum(num_steps + 1) * 2; 
            
            vec_Time      = NaT(1, n_approx, 'TimeZone', 'UTC');
            vec_SatID     = zeros(1, n_approx);
            vec_Range     = zeros(1, n_approx);
            vec_Elevation = zeros(1, n_approx);
            
            count = 0;
            
            for row = 1:height(intvls)
                sourceName = string(intvls.Source(row)); 
                sat_idx = sscanf(sourceName, "S4D_%d");
                
                t_start = intvls.StartTime(row);
                t_end   = intvls.EndTime(row);
                
                t_idx = find(simTimes >= t_start & simTimes <= t_end);
                
                if ~isempty(t_idx)
                    relevant_times = simTimes(t_idx);
                    num_new_points = numel(relevant_times);
                    
                    for t_k = 1:num_new_points
                        t_scalar = relevant_times(t_k);
                        [~, el, r] = aer(ue, sat(sat_idx), t_scalar);
                        
                        count = count + 1;
                        vec_Time(count)      = t_scalar;
                        vec_SatID(count)     = sat_idx;
                        vec_Range(count)     = r;
                        vec_Elevation(count) = el;
                    end
                end
            end
            
            if count > 0
                % Truncate to actual size
                raw_Time = vec_Time(1:count)';
                raw_SatID = vec_SatID(1:count)';
                raw_Range = vec_Range(1:count)';
                raw_Elev = vec_Elevation(1:count)';
                
                T = table(raw_Time, raw_SatID, raw_Range, raw_Elev, ...
                    'VariableNames', {'Time', 'SatID', 'Range', 'Elevation'});
                
                % Num satellites in view
                [G, ~] = findgroups(T.Time);
                num_visible = splitapply(@numel, T.SatID, G);
                
                T = sortrows(T, {'Time', 'Range'});
                
                % Save only the best (Closest) per timestamp
                [unique_times, idx_best] = unique(T.Time, 'first');
                T_best = T(idx_best, :);
                
                % Save filtered data
                SatData{pt} = struct(...
                    'Time',          T_best.Time', ...
                    'SatID',         T_best.SatID', ...
                    'Range',         T_best.Range', ...      
                    'Elevation_deg', T_best.Elevation', ...  
                    'Num_Visible',   num_visible', ...       
                    'UE_name',       ue_name_str, ...
                    'lat',           lat, ... % Note lowercase 'lat' to match link_calc expectation
                    'lon',           lon ...  % Note lowercase 'lon'
                );
            end
        end
        send(dq, []);
    end
    fprintf('\nSimulation Done. Saving results to %s...\n', dataFile);
    save(dataFile);
end
toc

%% Link Budget Configuration (From UE_focus_simulation)

% --- Downlink (Satellite -> UE) ---
Cfg.DL.Direction = "DL";
Cfg.DL.B     = 2e6;     % Downlink bandwidth
Cfg.DL.f     = 26e9;    % Downlink freq
Cfg.DL.P_tx  = 30;      % dBm (Satellite Power)
Cfg.DL.G_tx  = 36;      % dBi (Satellite Gain)
Cfg.DL.G_rx  = 36;      % dBi (UE Gain)
Cfg.DL.NF    = 5;       % dB  (UE Noise Figure)

% --- Uplink (UE -> Satellite) ---
Cfg.UL.Direction = "UL";
Cfg.UL.B     = 2e6;     % uplink bandwidth
Cfg.UL.f     = 28e9;    % uplink freq
Cfg.UL.P_tx  = 23;      % dBm (UE Power)
Cfg.UL.G_tx  = 36;      % dBi (UE Gain)
Cfg.UL.G_rx  = 36;      % dBi (Satellite Gain)
Cfg.UL.NF    = 3;       % dB  (Satellite Noise Figure)


%% Link Budget Calculations
fprintf('Starting Link Budget Calculations...\n');

% Note: 'SatData' here replaces 'UEs' array from the focus script
% iterating through all grid points (nPts)

% Preallocate fields to avoid structure inconsistency if needed, 
% but standard assignment works fine in Matlab for structs.

for pt = 1:nPts
    if ~isempty(SatData{pt})
        % Store SimData temporarily to match structure expected by link_calc
        % link_calc expects UE.SimData.Elevation_deg, etc.
        % Our SatData{pt} ALREADY has these fields directly. 
        % We will adapt link_calc slightly or construct a temp struct.
        
        % Construct a temp struct that mimics the 'UE' input for link_calc
        % link_calc accesses: UE.SimData.Elevation_deg, UE.lat, UE.lon
        
        % Let's act on the SatData{pt} directly, but we need to wrap it 
        % to make 'link_calc' work without modification, OR modify link_calc.
        % Modifying link_calc is cleaner given the context.
        % See 'modified_link_calc' below.
        
        SatData{pt}.DL = modified_link_calc(SatData{pt}, Cfg.DL);
        SatData{pt}.UL = modified_link_calc(SatData{pt}, Cfg.UL);
        
        % For backward compatibility with the plotting code that might expect 'Total_loss' 
        % directly on SatData{pt} or similar:
        % The plotting code in 'livescript_simulation' used:
        % SatData{pt}.Total_loss, SatData{pt}.SNR, etc.
        % We will copy the DL results to the top level or adjust the plotting code.
        % ADJUSTING PLOTTING CODE IS BETTER.
    end
end
fprintf('Link Budget Calculation complete.\n');


%% ---------------------------------------------------------
%  STATISTICS & VISUALIZATION MODULE
%  ---------------------------------------------------------

% 1. Setup & Selection
n_select = 5; 
selected_idxs = randperm(nPts, min(nPts, n_select));

% 2. Data Extraction
all_snr       = [];
all_thpt      = [];
all_loss      = [];
all_el        = [];
all_ranges_km = [];

fprintf('Extracting data for plots...\n');

for pt = 1:nPts
    if ~isempty(SatData{pt})
        % Using DOWNLINK (DL) for these statistics by default
        this_snr  = [SatData{pt}.DL.SNR];
        this_thpt = [SatData{pt}.DL.Throughput];
        this_loss = [SatData{pt}.DL.Total_loss];
        
        this_el   = [SatData{pt}.Elevation_deg];
        this_range= [SatData{pt}.Range];
        
        all_snr  = [all_snr, this_snr];
        all_thpt = [all_thpt, this_thpt];
        all_loss = [all_loss, this_loss];
        all_el   = [all_el, this_el];
        all_ranges_km = [all_ranges_km, this_range./1e3];
    end
end

%% --- FIGURE 1: INDIVIDUAL UE ANALYSIS (Random Sample) ---
figure('Name', 'Per-UE Performance (Random Sample)', 'Color', 'w');
t1 = tiledlayout(3,1, 'TileSpacing', 'compact');
title(t1, ['Per-UE Performance (DL)']);

% SNR
nexttile; hold on; grid on;
title('SNR Distribution');
xlabel('SNR (dB)'); ylabel('Count');
for idx = selected_idxs
    if ~isempty(SatData{idx})
        histogram(SatData{idx}.DL.SNR, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, ...
            'DisplayName', SatData{idx}.UE_name);
    end
end
legend('show', 'Location', 'eastoutside');

% Throughput
nexttile; hold on; grid on;
title('Throughput Distribution');
xlabel('Throughput (Mbps)'); ylabel('Count');
for idx = selected_idxs
    if ~isempty(SatData{idx})
        histogram(SatData{idx}.DL.Throughput./1e6, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, ...
            'DisplayName', SatData{idx}.UE_name);
    end
end

% Total Loss vs Elevation
nexttile; hold on; grid on;
title('Total Path Loss vs. Elevation');
xlabel('Elevation (deg)'); ylabel('Loss (dB)');
for idx = selected_idxs
    if ~isempty(SatData{idx})
        [sorted_el, sort_i] = sort(SatData{idx}.Elevation_deg);
        sorted_loss = SatData{idx}.DL.Total_loss(sort_i);
        plot(sorted_el, sorted_loss, '.-', 'LineWidth', 1, 'MarkerSize', 10, ...
            'DisplayName', SatData{idx}.UE_name);
    end
end
xlim([0 90]);


%% --- FIGURE 2: GLOBAL CONSTELLATION STATISTICS ---
figure('Name', 'Combined Constellation Stats', 'Color', 'w', 'Position', [100 100 1000 600]); 
t2 = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Global SNR Histogram
nexttile;
histogram(all_snr, 'Normalization', 'pdf', 'FaceColor', '#0072BD', 'EdgeColor', 'none');
grid on;
title('Global SNR PDF (DL)');
xlabel('SNR (dB)'); ylabel('Probability Density');

% Global Throughput Histogram
nexttile;
histogram(all_thpt./1e6, 'Normalization', 'pdf', 'FaceColor', '#D95319', 'EdgeColor', 'none');
grid on;
title('Global Throughput PDF (DL)');
xlabel('Throughput (Mbps)'); ylabel('Probability Density');

% Throughput CDF
nexttile;
[f_thpt, x_thpt] = ecdf(all_thpt./1e6);
plot(x_thpt, f_thpt, 'LineWidth', 2, 'Color', '#7E2F8E');
grid on;
title('Throughput CDF (DL)');
xlabel('Throughput (Mbps)'); ylabel('Probability  <=x');
xlim([0 max(x_thpt)]); 

% Summary Statistics Table
nexttile;
axis off; 
col1_str = {
    ['\bfSim Parameters (DL)\rm'];
    ['-------------------'];
    ['Total Datapoints: ' num2str(numel(all_snr))];
    ['Frequency:   ' num2str(Cfg.DL.f/1e9,'%.0f') ' GHz'];
    ['Bandwidth:   ' num2str(Cfg.DL.B/1e6,'%.2f') ' MHz'];
    ['Tx Power:    ' num2str(Cfg.DL.P_tx,'%.0f') ' dBm'];
};
col2_str = {
    ['\bfResults (DL)\rm'];
    ['-------------------'];
    ['\bfSNR (dB)\rm'];
    ['  Mean: ' num2str(mean(all_snr), '%.2f')];
    ['  Max:  ' num2str(max(all_snr), '%.2f')];
    ['\bfThroughput (Mbps)\rm'];
    ['  Mean: ' num2str(mean(all_thpt)/1e6, '%.2f')];
    ['  Max:  ' num2str(max(all_thpt)/1e6, '%.2f')];
};
text(0.05, 0.9, col1_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'Interpreter', 'tex');
text(0.55, 0.9, col2_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'Interpreter', 'tex');


%% --- LATITUDE ANALYSIS ---
% Setup Latitude Bins
lat_step = 5; 
lat_edges = -90:lat_step:90; 
lat_centers = lat_edges(1:end-1) + lat_step/2;
n_bins = length(lat_centers);

bin_counts = zeros(1, n_bins);
mean_Ar = zeros(1, n_bins); 
mean_Ag = zeros(1, n_bins); 
mean_Ac = zeros(1, n_bins); 
mean_As = zeros(1, n_bins); 
mean_At = zeros(1, n_bins); 

for pt = 1:nPts
    if ~isempty(SatData{pt})
        this_lat = SatData{pt}.lat;
        [~, bin_idx] = histc(this_lat, lat_edges);
        
        if bin_idx > 0 && bin_idx <= n_bins
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1;
            
            % Access DL Absorption struct
            abs_struct = SatData{pt}.DL.Absorption;
            
            mean_Ar(bin_idx) = mean_Ar(bin_idx) + max(abs_struct.Ar);
            mean_Ag(bin_idx) = mean_Ag(bin_idx) + max(abs_struct.Ag);
            mean_Ac(bin_idx) = mean_Ac(bin_idx) + max(abs_struct.Ac);
            mean_As(bin_idx) = mean_As(bin_idx) + max(abs_struct.As);
            mean_At(bin_idx) = mean_At(bin_idx) + max(abs_struct.At);
        end
    end
end

mask = bin_counts > 0;
mean_Ar(mask) = mean_Ar(mask) ./ bin_counts(mask);
mean_Ag(mask) = mean_Ag(mask) ./ bin_counts(mask);
mean_Ac(mask) = mean_Ac(mask) ./ bin_counts(mask);
mean_As(mask) = mean_As(mask) ./ bin_counts(mask);
mean_At(mask) = mean_At(mask) ./ bin_counts(mask);

valid_lats = lat_centers(mask);
valid_Ar = mean_Ar(mask);
valid_Ag = mean_Ag(mask);
valid_Ac = mean_Ac(mask);
valid_As = mean_As(mask);
valid_At = mean_At(mask);

figure('Name', 'Atmospheric Loss vs Latitude', 'Color', 'w');
axis; hold on; grid on;
Y_stack = [valid_Ag', valid_Ac', valid_As', valid_Ar'];
a = area(valid_lats, Y_stack);
a(1).FaceColor = '#A2142F'; a(1).FaceAlpha = 0.4; 
a(2).FaceColor = '#77AC30'; a(2).FaceAlpha = 0.4; 
a(3).FaceColor = '#EDB120'; a(3).FaceAlpha = 0.4; 
a(4).FaceColor = '#0072BD'; a(4).FaceAlpha = 0.4; 
p_total = plot(valid_lats, valid_At, 'k-o', 'LineWidth', 2.5, 'MarkerFaceColor', 'k', 'MarkerSize', 4);
title({'Atmospheric Loss Contributions by Latitude (DL)', '99% Availability'});
xlabel('Latitude (deg)'); ylabel('Max Path Loss (dB)');
legend([a(4), a(3), a(2), a(1), p_total], {'Rain', 'Scintillation', 'Clouds', 'Gases', 'Total Model'}, 'Location', 'best');

%% --- MAP VISUALIZATIONS ---
% Compute stats for map plotting
minNumberSatellites = zeros(nPts, 1);
meanNumberSatellites = zeros(nPts, 1);
lat_vector = zeros(nPts, 1);
lon_vector = zeros(nPts, 1);

for pt = 1:nPts
    if ~isempty(SatData{pt})
        minNumberSatellites(pt) = min(SatData{pt}.Num_Visible);
        meanNumberSatellites(pt) = mean(SatData{pt}.Num_Visible);
    else
        minNumberSatellites(pt) = 0;
        meanNumberSatellites(pt) = 0;
    end
    lat_vector(pt) = LatGrid(pt);
    lon_vector(pt) = LonGrid(pt);
end

% Example Map: Mean Number of Satellites
figure;
ax = axesm('lambertstd', 'MapLatLimit', [min(lat_vec) max(lat_vec)], 'MapLonLimit', [min(lon_vec) max(lon_vec)], 'Frame', 'on', 'Grid', 'on');
axis off;
land = shaperead('landareas.shp','UseGeoCoords',true);
geoshow([land.Lat], [land.Lon], 'DisplayType','polygon', 'FaceColor',[0.8 0.8 0.8]);

% Interpolation
lat_lim = [min(lat_vec) max(lat_vec)];
lon_lim = [min(lon_vec) max(lon_vec)];
[LonG, LatG] = meshgrid(linspace(lon_lim(1), lon_lim(2), 200), linspace(lat_lim(1), lat_lim(2), 200));
ValG = griddata(lon_vector, lat_vector, meanNumberSatellites, LonG, LatG, 'cubic');

surfm(LatG, LonG, ValG,'FaceAlpha',0.5);
cb = colorbar; ylabel(cb, 'Mean Number of Satellites');
title('Mean Number of Satellites Coverage');


%% --- HELPER FUNCTIONS ---

function Link = modified_link_calc(UE_struct, cfg)
    % Adapted to take a plain struct with vectors
    
    el_vec    = [UE_struct.Elevation_deg];
    range_vec = [UE_struct.Range];
    lat       = UE_struct.lat;
    lon       = UE_struct.lon;

    Link.Frequency = cfg.f;
    Link.Bandwidth = cfg.B;

    % losses
    Link.FSPL = FSPL_calc(cfg.f, range_vec);
    [Link.Absorption, Link.T_antenna] = Absorption_calc(cfg.f, el_vec, lat, lon);
    Link.Steering_loss = Steering_loss_calc(el_vec, range_vec);
    Link.Total_loss = Link.FSPL + Link.Absorption.At + Link.Steering_loss;
    
    if (cfg.Direction == "UL")
        Link.T_antenna = zeros(size(Link.T_antenna)) + 290; % Overwrite antenna temp by 290K for UL
    end

    % power and noise
    Link.P_noise = Noise_density_calc(Link.T_antenna, cfg.NF) + 10*log10(cfg.B);
    Link.Rx_Power = cfg.P_tx + cfg.G_tx + cfg.G_rx - Link.Total_loss;

    % link quality
    Link.SNR = Link.Rx_Power - Link.P_noise;
    Link.Throughput = Throughput_calc(Link.SNR, cfg.B);
end

function loss_dB = FSPL_calc(f, range_vec)
    c = physconst('LightSpeed');
    lambda = c / f;
    loss_dB = 10 * log10(((4 * pi * range_vec) / lambda).^2);
end

function loss_dB = Steering_loss_calc(el_vec, range_vec)
    Re = 6378.14e3;     % Earth Radius
    cos_exponent = 1.5; % Default scan loss parameter

    theta_ue = 90 - el_vec; 
    
    R_sat = sqrt(Re^2 + range_vec.^2 + 2 * Re .* range_vec .* sind(el_vec));
    sin_theta_sat = (Re .* cosd(el_vec)) ./ R_sat;
    theta_sat = asind(sin_theta_sat);
    
    loss_factor = (cosd(theta_ue).^cos_exponent) .* ...
                  (cosd(theta_sat).^cos_exponent);

    loss_dB = -10 * log10(loss_factor); 
end

function N0_dBmHz = Noise_density_calc(antenna_temp, NF)
    k_dB = 10*log10(1.38e-23); 
    F_lin = 10^(NF/10);
    T_rx = (F_lin - 1) * 290;
    T_sys = antenna_temp + T_rx;
    N0_dBW_Hz = 10*log10(T_sys) + k_dB;
    N0_dBmHz = N0_dBW_Hz + 30;
end

function Throughput = Throughput_calc(SNR_dB,B)
    BW_eff = 0.56;
    eta = 1;
    SNR_eff_dB = 2;

    SNR_eff = 10^(SNR_eff_dB/10);
    SNR = 10.^(SNR_dB/10);

    Throughput = B*BW_eff*eta*log2(1+SNR./SNR_eff);
end

function [loss_struct, sky_temp_K] = Absorption_calc(f, el_vec, lat, lon)
    el_grid = 20:10:90; 
    n_grid = length(el_grid);
    
    grid_Ag = zeros(1, n_grid);
    grid_Ac = zeros(1, n_grid);
    grid_Ar = zeros(1, n_grid);
    grid_As = zeros(1, n_grid);
    grid_At = zeros(1, n_grid);
    grid_Tsky = zeros(1, n_grid);
    
    prob = 1; 
    pol = 45; 
    Antenna_Dia = 0.5; 
    Antenna_eff = 0.5; 
    
    for i = 1:n_grid
        cfg = p618Config(...
            'Frequency', f, ...
            'ElevationAngle', el_grid(i), ... 
            'Latitude', lat, ...
            'Longitude', lon, ...
            'TotalAnnualExceedance', prob, ... 
            'PolarizationTiltAngle', pol, ...
            'AntennaDiameter', Antenna_Dia, ...     
            'AntennaEfficiency', Antenna_eff);      
        
        [pl, ~, tsky] = p618PropagationLosses(cfg);
        
        grid_Ag(i) = pl.Ag;
        grid_Ac(i) = pl.Ac;
        grid_Ar(i) = pl.Ar;
        grid_As(i) = pl.As;
        grid_At(i) = pl.At;
        grid_Tsky(i) = tsky;
    end

    loss_struct.Ag = interp1(el_grid, grid_Ag, el_vec, 'linear');
    loss_struct.Ac = interp1(el_grid, grid_Ac, el_vec, 'linear');
    loss_struct.Ar = interp1(el_grid, grid_Ar, el_vec, 'linear');
    loss_struct.As = interp1(el_grid, grid_As, el_vec, 'linear');
    loss_struct.At = interp1(el_grid, grid_At, el_vec, 'linear');
    sky_temp_K = interp1(el_grid, grid_Tsky, el_vec, 'linear');
end

function updateCommandLine(total_pts)
    persistent p
    if isempty(p) || p == total_pts
        p = 0;
    end
    p = p + 1;
    percent = (p / total_pts) * 100;
    if p > 1
        fprintf(repmat('\b', 1, 7)); 
    end
    fprintf('%6.2f%%', percent);
end
