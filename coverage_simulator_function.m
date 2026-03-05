function metrics = coverage_simulator_function(Cfg,plot_results, use_parallel)
% RUN_SATELLITE_SIM Simulates satellite coverage and link budget.
% Inputs: Cfg - Struct containing all configuration parameters
% Outputs: metrics - Struct containing key performance indicators for optimization

    % fprintf('\n======================================================\n');
    fprintf('\n Starting Simulation: %d Sats, %.1f deg Inclination\n', Cfg.Total_sats, Cfg.Inclination);
    % fprintf('======================================================\n');

    if plot_results
        %% Create Output Directory
        % Format: simulation_output/orbit_sats_inclination_phasing_date
        date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name = sprintf('%.0f_%d_%.0f_%.1f_%s', ...
            Cfg.Orbit_height/1e3, Cfg.Total_sats, Cfg.Inclination, Cfg.Phasing, date_str);
        out_dir = fullfile('simulation_output', folder_name);
        
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        % fprintf('Results will be saved to: %s\n', out_dir);
    
        % Save the configuration file immediately
        save(fullfile(out_dir, 'Config.mat'), 'Cfg');
    end

    %% Scenario & Constellation Setup
    % fprintf('Defining constellation...\n');
    sc = satelliteScenario;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;
    
    simTimes = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
    simTimes.TimeZone = 'UTC';

    r_earth = 6378.14e3;
    if Cfg.WalkerStar == true
        % CUSTOM ASYMMETRICAL WALKER STAR
        a = r_earth + Cfg.Orbit_height;
        e = 0; % Circular orbit
        inc = Cfg.Inclination;
        argPer = 0; 
        
        S = Cfg.Total_sats / Cfg.Num_planes; % Sats per plane
        P = Cfg.Num_planes; % Keep as total planes over 180 degrees
        

        corotating_gap = (180 - Cfg.Seam_gap) / (P - 1);
        
        % Calculate the Walker Phasing offset
        phase_shift = (Cfg.Phasing * 360) / Cfg.Total_sats;
        
        sat_array = []; % Temporary array to hold satellites
        
        for p = 1:P
            % Calculate asymmetrical RAAN
            raan = (p - 1) * corotating_gap;
            
            for s = 1:S
                % Calculate True Anomaly (position in the orbit + staggering)
                nu = (s - 1) * (360 / S) + (p - 1) * phase_shift;
                nu = mod(nu, 360); % Keep it cleanly within 0-360 degrees
                
                % Must match your existing "S4D_X" naming convention
                sat_num = (p - 1) * S + s;
                name = sprintf("S4D_%d", sat_num);
                
                % Create satellite and add to scenario
                new_sat = satellite(sc, a, e, inc, raan, argPer, nu, ...
                    "Name", name, "OrbitPropagator", "sgp4");
                
                sat_array = [sat_array, new_sat]; 
            end
        end
        sat = sat_array; % Assign array to the expected variable
    else
        sat = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, ...
        Cfg.Total_sats, ...
        Cfg.Num_planes, ...
        Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");
    end

    %% UE Grid Setup
    [LonGrid, LatGrid] = meshgrid(Cfg.Lon_vec, Cfg.Lat_vec);
    NumUEs = numel(LonGrid);
    UEs = cell(NumUEs, 1);
    for idx = 1:NumUEs
        UEs{idx}.Lat = LatGrid(idx);
        UEs{idx}.Lon = LonGrid(idx);
        UEs{idx}.Name = sprintf('UE%dLat%.0f', idx, UEs{idx}.Lat);
    end

    %% Coverage Simulation
    % fprintf('Defining UEs and calculating links...\n');
    
    dq = parallel.pool.DataQueue;
    % Call it once with the reset flag to print the 0.00%
    updateLiveScriptProgress(NumUEs, true); 
    % Tell the queue to run it without the flag after each loop
    afterEach(dq, @(~) updateLiveScriptProgress(NumUEs, false));

    if use_parallel
        num_workers = Inf; % Use the active parallel pool
    else
        num_workers = 0;   % Force serial execution (acts exactly like a 'for' loop)
    end


    tic
    parfor (idx = 1:NumUEs,num_workers)
        current_UE = UEs{idx};
        ue = groundStation(sc, current_UE.Lat, current_UE.Lon, ...
            'Name', current_UE.Name, 'MinElevationAngle', Cfg.Min_elevation_UE);
        
        ac = access(sat, ue);
        intvls = accessIntervals(ac);
        
        [~, el_mat, r_mat, sim_Times] = aer(ue, sat);
        sim_Times = sim_Times'; 
        nT = length(sim_Times);
        
        valid_mask = false(size(el_mat));
        
        for row = 1:height(intvls)
            sourceName = string(intvls.Source(row));
            sat_idx = sscanf(sourceName, "S4D_%d");
            
            delta_start = seconds(intvls.StartTime(row) - sc.StartTime);
            delta_end   = seconds(intvls.EndTime(row)   - sc.StartTime);
            
            idx_start = round(delta_start / sc.SampleTime) + 1;
            idx_end   = round(delta_end   / sc.SampleTime) + 1;
            valid_mask(sat_idx, idx_start:idx_end) = true;
        end
        
        Num_visible = sum(valid_mask,1); 
        has_service = Num_visible > 0;
        
        r_temp = r_mat;
        r_temp(~valid_mask) = Inf; 
        [best_ranges, best_sat_idx] = min(r_temp, [], 1); 
        
        final_Range = nan(1, nT);
        final_El    = nan(1, nT);
        final_SatID = nan(1, nT);
        
        if any(has_service)
            final_Range(has_service) = best_ranges(has_service);
            best_sats_valid = best_sat_idx(has_service);
            final_SatID(has_service) = best_sats_valid;
            
            valid_cols = find(has_service);
            num_rows = size(el_mat, 1);
            lin_idxs = best_sats_valid + (valid_cols - 1) * num_rows;
            final_El(has_service) = el_mat(lin_idxs); % linear index needed for a 1d output
        end
        
        current_UE.SimData = struct(... 
            'Time',          sim_Times', ...            
            'SatID',         final_SatID, ...          
            'Range',         final_Range, ...          
            'Elevation_deg', final_El, ...             
            'Num_visible',   Num_visible ...           
        );
        
        durations = seconds(intvls.EndTime - intvls.StartTime);
        current_UE.SimData.TimeStats = struct(...
            'Mean',   mean(durations), ...
            'Max',    max(durations), ...
            'Min',    min(durations), ...
            'Total',  sum(durations), ...
            'Num_Passes', numel(durations) ... 
        );
        
        UEs{idx} = current_UE;
        send(dq, []);
    end
    fprintf('\nGeometry calculation complete (%.1f sec).', toc);

    %% Link Budget Calculation

    updateLiveScriptProgress(NumUEs, true); 
    % Tell the queue to run it without the flag after each loop
    afterEach(dq, @(~) updateLiveScriptProgress(NumUEs, false));
    tic
    parfor (idx = 1:NumUEs,num_workers)
        UEs{idx}.DL = link_calc(UEs{idx}, Cfg.DL);
        send(dq, []);
    end
    fprintf('\n Loss calculation complete (%.1f sec).', toc);

    %% Data Aggregation & Metrics
    all_snr  = [];
    all_thpt = [];
    all_loss = [];
    nT = length(UEs{1}.SimData.Time);
    counts = zeros(NumUEs, nT);
    
    for idx = 1:NumUEs
        all_snr  = [all_snr, [UEs{idx}.DL.SNR]];
        all_thpt = [all_thpt, [UEs{idx}.DL.Throughput]];
        all_loss = [all_loss, [UEs{idx}.DL.Total_loss]];
        counts(idx,:) = UEs{idx}.SimData.Num_visible;
    end

    % Coverage Stats
    more_1_satellites = (counts >= 1);
    prob_coverage = 100 * sum(more_1_satellites, 2) ./ nT; % Array per UE
    minNumberSatellites = min(counts, [], 2);
    meanNumberSatellites = mean(counts, 2);
    
    maxGapMinutes = zeros(NumUEs,1);
    for idx = 1:NumUEs
        row = (counts(idx,:) == 0);  
        maxZeroStreak = 0; currentStreak = 0;
        for k = 1:length(row)
            if row(k)
                currentStreak = currentStreak + 1;
                maxZeroStreak = max(maxZeroStreak, currentStreak);
            else
                currentStreak = 0;
            end
        end
        maxGapMinutes(idx) = maxZeroStreak * Cfg.SampleTime / 60;
    end

    % Compile specific return metrics for the search algorithm
    metrics.worst_coverage_percent = min(prob_coverage);
    metrics.worst_gap_minutes      = max(maxGapMinutes);
    metrics.throughput_10pct_Mbps  = prctile(all_thpt, 10) / 1e6;
    metrics.throughput_mean_Mbps   = mean(all_thpt, 'omitnan') / 1e6;
    metrics.Num_visible = counts; % returns num visible timeseries

    fprintf('Metrics Calculated: Min Cov = %.2f%%, Max Gap = %.1f min, 10%% Thpt = %.2f Mbps\n', ...
        metrics.worst_coverage_percent, metrics.worst_gap_minutes, metrics.throughput_10pct_Mbps);

    %% Plot Generation & Saving
    % fprintf('Generating and saving plots...\n');
    if plot_results
        
        % Combined Constellation Stats
        f1 = figure('Visible', 'off', 'Name', 'Combined Constellation Stats', 'Color', 'w', 'Position', [100 100 1000 600]); 
        t2 = tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Tile 1: Throughput PDF
        nexttile; histogram(all_thpt./1e6, 'Normalization', 'pdf', 'FaceColor', '#D95319', 'EdgeColor', 'none');
        grid on; title('Global Throughput PDF (DL)'); xlabel('Throughput (Mbps)'); ylabel('Probability Density');
        
        % Tile 2: Throughput CDF
        nexttile; [f_thpt, x_thpt] = ecdf(all_thpt./1e6); plot(x_thpt, f_thpt, 'LineWidth', 2, 'Color', '#7E2F8E');
        grid on; title('Throughput CDF (DL)'); xlabel('Throughput (Mbps)'); ylabel('Probability \leq x'); xlim([0 max(x_thpt)]);
        
        % Tile 3 & 4 (Merged to span the whole bottom row for our text)
        nexttile(3, [1 2]); axis off; 
        
        % Determine Walker string
        if Cfg.WalkerStar
            walker_str = 'Walker Star';
        else
            walker_str = 'Walker Delta';
        end

        % Column 1: Constellation Data
        col1_str = {
            '\bfConstellation Settings\rm';
            ['Type:         ' walker_str];
            ['Total Sats:   ' num2str(Cfg.Total_sats)];
            ['Planes/Sats:  ' num2str(Cfg.Num_planes) ' / ' num2str(Cfg.Sats_per_plane)];
            ['Inclination:  ' num2str(Cfg.Inclination, '%.1f') '\circ'];
            ['Phasing:      ' num2str(Cfg.Phasing)];
            ['Altitude:     ' num2str(Cfg.Orbit_height/1e3, '%.0f') ' km'];
            ['Min Elev:     ' num2str(Cfg.Min_elevation_UE, '%.1f') '\circ'];
        };
        
        % Column 2: Link Budget / RF Data
        col2_str = {
            '\bfLink Budget Specs\rm';
            ['Direction:    ' char(Cfg.DL.Direction)];
            ['Freq / BW:    ' num2str(Cfg.DL.f/1e9, '%.2f') ' GHz / ' num2str(Cfg.DL.B/1e6, '%.1f') ' MHz'];
            ['Tx Type/Gain: ' char(Cfg.DL.Tx_type) ' / ' num2str(Cfg.DL.G_tx, '%.1f') ' dBi'];
            ['P_tx / EIRP:  ' num2str(Cfg.DL.P_tx, '%.1f') ' W / ' num2str(Cfg.DL.EIRP_dBm, '%.1f') ' dBm'];
            ['Rx Type/Gain: ' char(Cfg.DL.Rx_type) ' / ' num2str(Cfg.DL.G_rx, '%.1f') ' dBi'];
            ['Noise Fig:    ' num2str(Cfg.DL.NF, '%.1f') ' dB'];
        };
        
        % Column 3: Performance Results
        col3_str = {
            '\bfPerformance Results\rm';
            ['10% Throughput:  ' num2str(metrics.throughput_10pct_Mbps, '%.2f') ' Mbps'];
            ['Mean Throughput: ' num2str(metrics.throughput_mean_Mbps, '%.2f') ' Mbps'];
            ['Worst Coverage:  ' num2str(metrics.worst_coverage_percent, '%.2f') ' %'];
            ['\bfSNR (dB)\rm'];
            ['  Mean: ' num2str(mean(all_snr), '%.2f')];
            ['  Min:  ' num2str(min(all_snr), '%.2f')];
            ['  Max:  ' num2str(max(all_snr), '%.2f')];
        };
        
        % Plot text in 3 evenly spaced columns
        text(0.05, 0.95, col1_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        text(0.38, 0.95, col2_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        text(0.72, 0.95, col3_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        
        exportgraphics(f1, fullfile(out_dir, 'Global_Stats.png'), 'Resolution', 300);
        close(f1);

        % 
        % 
        % 
        % 
        % 
        % % Combined Constellation Stats
        % f1 = figure('Visible', 'off', 'Name', 'Combined Constellation Stats', 'Color', 'w', 'Position', [100 100 1000 600]); 
        % t2 = tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        % 
        % nexttile; histogram(all_snr, 'Normalization', 'pdf', 'FaceColor', '#0072BD', 'EdgeColor', 'none');
        % grid on; title('Global SNR PDF (DL)'); xlabel('SNR (dB)'); ylabel('Probability Density');
        % 
        % nexttile; histogram(all_thpt./1e6, 'Normalization', 'pdf', 'FaceColor', '#D95319', 'EdgeColor', 'none');
        % grid on; title('Global Throughput PDF (DL)'); xlabel('Throughput (Mbps)'); ylabel('Probability Density');
        % 
        % nexttile; [f_thpt, x_thpt] = ecdf(all_thpt./1e6); plot(x_thpt, f_thpt, 'LineWidth', 2, 'Color', '#7E2F8E');
        % grid on; title('Throughput CDF (DL)'); xlabel('Throughput (Mbps)'); ylabel('Probability <= x'); xlim([0 max(x_thpt)]);
        % 
        % nexttile; axis off; 
        % col1_str = {
        %     ['\bfConstellation\rm'];
        %     ['Total Satellites: ' num2str(Cfg.Total_sats)];
        %     ['Inclination:      ' num2str(Cfg.Inclination, '%.1f') '\circ'];
        %     ['Altitude:         ' num2str(Cfg.Orbit_height/1e3, '%.0f') ' km'];
        % };
        % col2_str = {
        %     ['\bfDL Performance Results\rm'];
        %     ['10% Throughput: ' num2str(metrics.throughput_10pct_Mbps, '%.2f') ' Mbps'];
        %     ['Mean Throughput: ' num2str(metrics.throughput_mean_Mbps, '%.2f') ' Mbps'];
        %     ['Worst Coverage: ' num2str(metrics.worst_coverage_percent, '%.2f') ' %'];
        % };
        % text(0.05, 0.9, col1_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        % text(0.55, 0.9, col2_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        % 
        % exportgraphics(f1, fullfile(out_dir, 'Global_Stats.png'), 'Resolution', 300);
        % close(f1);
    
        % Mapping (Probability of Service)
        lat_vector = cellfun(@(x) x.Lat, UEs)';
        lon_vector = cellfun(@(x) x.Lon, UEs)';
        
        % 1) Geographical Grid (Calculated once for all maps)
        lat_lim = [min(Cfg.Lat_vec) max(Cfg.Lat_vec)];
        lon_lim = [min(Cfg.Lon_vec) max(Cfg.Lon_vec)];
        nLat = 200; nLon = 200;
        [LonG, LatG] = meshgrid(linspace(lon_lim(1), lon_lim(2), nLon), linspace(lat_lim(1), lat_lim(2), nLat));
        
        try
            %% MAP 1: Min Number of Satellites
            ValG_min = griddata(lon_vector, lat_vector, minNumberSatellites, LonG, LatG, 'cubic');
            
            f2 = figure('Visible', 'off', 'Color', 'w');
            ax2 = axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, ...
                        'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
            axis off;  
            land = shaperead('landareas.shp', 'UseGeoCoords', true);
            geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
            surfm(LatG, LonG, ValG_min, 'FaceAlpha', 0.5);
            
            for i = 1:NumUEs
                textm(lat_vector(i), lon_vector(i), sprintf('%d', minNumberSatellites(i)), ...
                      'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                      'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k'); 
            end
            
            cb2 = colorbar; caxis([0 max([minNumberSatellites; 1])]); ylabel(cb2, 'Min Number of Satellites');
            set(gca, 'FontSize', 20);
            exportgraphics(f2, fullfile(out_dir, 'Map_Min_Sats.png'), 'Resolution', 300);
            close(f2);
            
            %% MAP 2: Mean Number of Satellites
            ValG_mean = griddata(lon_vector, lat_vector, meanNumberSatellites, LonG, LatG, 'cubic');
            
            f3 = figure('Visible', 'off', 'Color', 'w');
            ax3 = axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, ...
                        'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
            axis off;  
            geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
            surfm(LatG, LonG, ValG_mean, 'FaceAlpha', 0.5);
            
            for i = 1:NumUEs
                textm(lat_vector(i), lon_vector(i), sprintf('%.1f', meanNumberSatellites(i)), ...
                      'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                      'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
            end
            
            cb3 = colorbar; caxis([0 max([meanNumberSatellites; 1])]); ylabel(cb3, 'Mean Number of Satellites');
            set(gca, 'FontSize', 20);
            exportgraphics(f3, fullfile(out_dir, 'Map_Mean_Sats.png'), 'Resolution', 300);
            close(f3);
            
            %% MAP 3: Probability of Coverage (>= 1 Satellite)
            ValG_prob = griddata(lon_vector, lat_vector, prob_coverage, LonG, LatG, 'cubic');
            
            f4 = figure('Visible', 'off', 'Color', 'w');
            ax4 = axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, ...
                        'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
            axis off; 
            geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
            surfm(LatG, LonG, ValG_prob, 'FaceAlpha', 0.5);
            
            for i = 1:NumUEs
                textm(lat_vector(i), lon_vector(i), sprintf('%.0f', prob_coverage(i)), ...
                      'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                      'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k'); 
            end
            
            cb4 = colorbar; caxis([0 100]); ylabel(cb4, 'P(Sat \geq 1) [%]');
            set(gca, 'FontSize', 20);
            exportgraphics(f4, fullfile(out_dir, 'Map_Coverage_Prob.png'), 'Resolution', 300);
            close(f4);
            
        catch ME
            fprintf('Warning: Map generation failed. Ensure Mapping Toolbox is installed and grid is not too sparse.\n');
            disp(ME.message);
        end
    
        fprintf('Run Complete. All data saved to %s\n', out_dir);
    end
end

%% ========================================================================
% HELPER FUNCTIONS (Kept internally so you only need one file)
% =========================================================================

function Link = link_calc(UE, cfg)
    el_vec    = [UE.SimData.Elevation_deg];
    range_vec = [UE.SimData.Range];
    lat       = UE.Lat;
    lon       = UE.Lon;
    Link.Frequency = cfg.f;
    Link.Bandwidth = cfg.B;
    
    Link.FSPL = FSPL_calc(cfg.f, range_vec);
    [Link.Absorption, Link.T_antenna] = Absorption_calc(cfg.f, el_vec, lat, lon);
    Link.Steering_loss = Steering_loss_calc(el_vec, range_vec, cfg.Tx_type, cfg.Rx_type);
    Link.Total_loss = Link.FSPL + Link.Absorption.At + Link.Steering_loss;
    
    if (isfield(cfg, 'Direction') && cfg.Direction == "UL")
        Link.T_antenna = zeros(size(Link.T_antenna)) + 290;
    end
    
    Link.P_noise = Noise_density_calc(Link.T_antenna, cfg.NF) + 10*log10(cfg.B);
    Link.Rx_Power = cfg.EIRP_dBm + cfg.G_rx - Link.Total_loss;
    Link.SNR = Link.Rx_Power - Link.P_noise;
    
    Link.Throughput = nan(size(Link.SNR));
    valid_idx = ~isnan(Link.SNR);
    Link.Throughput(valid_idx) = Throughput_calc(Link.SNR(valid_idx), cfg.B);
    Link.Throughput(~valid_idx) = 0; % Force 0 for NaN
end

function loss_dB = FSPL_calc(f, range_vec)
    c = physconst('LightSpeed');
    lambda = c / f;
    loss_dB = 10 * log10(((4 * pi * range_vec) / lambda).^2);
end

function Throughput = Throughput_calc(SNR_dB,B)
    BW_eff = 0.56;
    eta = 1;
    SNR_eff_dB = 2;
    SNR_eff = 10^(SNR_eff_dB/10);
    SNR = 10.^(SNR_dB/10);
    Throughput = B*BW_eff*eta.*log2(1+SNR./SNR_eff);
end

function [loss_struct, sky_temp_K] = Absorption_calc(f, el_vec, lat, lon)
    el_grid = 20:10:90; 
    n_grid = length(el_grid);
    grid_Ag = zeros(1, n_grid); grid_Ac = zeros(1, n_grid);
    grid_Ar = zeros(1, n_grid); grid_As = zeros(1, n_grid);
    grid_At = zeros(1, n_grid); grid_Tsky = zeros(1, n_grid);
    
    for i = 1:n_grid
        cfg = p618Config('Frequency', f, 'ElevationAngle', el_grid(i), ... 
            'Latitude', lat, 'Longitude', lon, 'TotalAnnualExceedance', 1, ... 
            'PolarizationTiltAngle', 45, 'AntennaDiameter', 0.5, 'AntennaEfficiency', 0.5);      
        [pl, ~, tsky] = p618PropagationLosses(cfg);
        grid_Ag(i) = pl.Ag; grid_Ac(i) = pl.Ac; grid_Ar(i) = pl.Ar;
        grid_As(i) = pl.As; grid_At(i) = pl.At; grid_Tsky(i) = tsky;
    end
    
    loss_struct.Ag = interp1(el_grid, grid_Ag, el_vec, 'linear');
    loss_struct.Ac = interp1(el_grid, grid_Ac, el_vec, 'linear');
    loss_struct.Ar = interp1(el_grid, grid_Ar, el_vec, 'linear');
    loss_struct.As = interp1(el_grid, grid_As, el_vec, 'linear');
    loss_struct.At = interp1(el_grid, grid_At, el_vec, 'linear');
    sky_temp_K = interp1(el_grid, grid_Tsky, el_vec, 'linear');
end

function loss_dB = Steering_loss_calc(el_vec, range_vec, TxType, RxType)
    Re = 6378.14e3;     
    cos_exponent = 1.5; 
    theta_rx = 90 - el_vec; 
    R_sat_sq = Re^2 + range_vec.^2 + 2 * Re .* range_vec .* sind(el_vec);
    R_sat = sqrt(R_sat_sq);
    sin_theta_tx = (Re .* cosd(el_vec)) ./ R_sat;
    theta_tx = asind(sin_theta_tx); 
    
    loss_tx = 0; loss_rx = 0;
    if contains(TxType, 'array')
        loss_tx = -10 * log10(cosd(theta_tx).^cos_exponent);
    end
    if contains(RxType, 'array')
        loss_rx = -10 * log10(cosd(theta_rx).^cos_exponent);
    end
    loss_dB = loss_tx + loss_rx;
end

function N0_dBmHz = Noise_density_calc(antenna_temp, NF)
    k_dB = 10*log10(1.38e-23); 
    F_lin = 10^(NF/10);
    T_rx = (F_lin - 1) * 290;
    T_sys = antenna_temp + T_rx;
    N0_dBmHz = (10*log10(T_sys) + k_dB) + 30;
end

function updateLiveScriptProgress(total_pts, reset_flag)
    persistent p last_percent
    
    % Initialization / Reset
    if nargin > 1 && reset_flag
        p = 0;
        last_percent = 0;
        % fprintf('Progress: 0%%... ');
        return;
    end
    
    if isempty(p)
        p = 0;
        last_percent = 0;
    end
    
    p = p + 1;
    current_percent = floor((p / total_pts) * 100);
    
    % Only print to the Live Script when we cross a 10% threshold
    if current_percent >= last_percent + 10
        fprintf('%d%%... ', current_percent);
        last_percent = current_percent;
    end
    
    % Cap it off cleanly when finished
    if p >= total_pts
        fprintf('Done!\n');
        p = 0; % Auto-reset for the next run
    end
end