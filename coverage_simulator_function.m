function metrics = coverage_simulator_function(Cfg,plot_results, use_parallel, calc_link)
% RUN_SATELLITE_SIM Simulates satellite coverage and link budget.
% Inputs: Cfg - Struct containing all configuration parameters
% Outputs: metrics - Struct containing key performance indicators for optimization

    % fprintf('\n======================================================\n');
    fprintf('\n Starting Simulation: %d Sats, %.1f deg Inclination\n', Cfg.Total_sats, Cfg.Inclination);
    % fprintf('======================================================\n');
    

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
        sats = asymmetrical_walker_star_generation(Cfg.Orbit_height, Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane);
    else
        sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, ...
        Cfg.Total_sats, ...
        Cfg.Num_planes, ...
        Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");
    end

    %% Create the UEs Array
    if Cfg.Equal_UE_area == true
        [UE_lats, UE_lons] = generate_equal_ish_area_UEs(Cfg.Lat_vec, Cfg.Lon_vec);
    elseif Cfg.Accept_Flat_UE_array == true
        UE_lats = Cfg.Flat_UE_array.Lats;
        UE_lons = Cfg.Flat_UE_array.Lons;
    else
        [UE_lats, UE_lons] = meshgrid(Cfg.Lat_vec, Cfg.Lon_vec);
    end
    NumUEs = length(UE_lats);
    UEs = cell(NumUEs, 1);
    
    for idx = 1:NumUEs
        UEs{idx}.Lat = UE_lats(idx);
        UEs{idx}.Lon = UE_lons(idx);
        UEs{idx}.Name = sprintf('UE%d', idx);
    end

    %% Coverage Simulation
    % fprintf('Defining UEs and calculating links...\n');
    
    dq = parallel.pool.DataQueue;
    updateLiveScriptProgress(NumUEs, true); 
    afterEach(dq, @(~) updateLiveScriptProgress(NumUEs, false));

    if use_parallel
        num_workers = Inf; % Use the active parallel pool
    else
        num_workers = 0;   % acts exactly like a 'for' loop
    end

    min_elevation_UE = Cfg.Min_elevation_UE;
    tic
    parfor (idx = 1:NumUEs,num_workers)
        current_UE = UEs{idx};
        ue = groundStation(sc, current_UE.Lat, current_UE.Lon, ...
            'Name', current_UE.Name, 'MinElevationAngle', min_elevation_UE);
        
        
        
        [az_mat, el_mat, r_mat, sim_Times] = aer(ue, sats);
        sim_Times = sim_Times'; 
        nT = length(sim_Times);
        
        valid_mask = el_mat >= min_elevation_UE;


        % ac = access(sats, ue);
        % intvls = accessIntervals(ac);
        % valid_mask = false(size(el_mat));
        % 
        % for row = 1:height(intvls)
        %     sourceName = string(intvls.Source(row));
        %     sat_idx = sscanf(sourceName, "S4D_%d");
        % 
        %     delta_start = seconds(intvls.StartTime(row) - sc.StartTime);
        %     delta_end   = seconds(intvls.EndTime(row)   - sc.StartTime);
        % 
        %     idx_start = round(delta_start / sc.SampleTime) + 1;
        %     idx_end   = round(delta_end   / sc.SampleTime) + 1;
        %     valid_mask(sat_idx, idx_start:idx_end) = true;
        % end

        Num_visible = sum(valid_mask,1); 
        has_service = Num_visible > 0;

        r_temp = r_mat;
        r_temp(~valid_mask) = Inf; 
        [best_ranges, best_sat_idx] = min(r_temp, [], 1); 

        final_Range = nan(1, nT);
        final_El    = nan(1, nT);
        final_Az    = nan(1, nT);
        final_SatID = nan(1, nT);

        if any(has_service)
            final_Range(has_service) = best_ranges(has_service);
            best_sats_valid = best_sat_idx(has_service);
            final_SatID(has_service) = best_sats_valid;

            valid_cols = find(has_service);
            num_rows = size(el_mat, 1);
            lin_idxs = best_sats_valid + (valid_cols - 1) * num_rows;
            final_El(has_service) = el_mat(lin_idxs); % linear index needed for a 1d output
            final_Az(has_service) = az_mat(lin_idxs);
        end
        
        current_UE.SimData = struct(... 
            'Time',          sim_Times', ...            
            'SatID',         final_SatID, ...          
            'Range',         final_Range, ...          
            'Elevation_deg', final_El, ...
            'Azimuth_deg',   final_Az, ...
            'Num_visible',   Num_visible ...           
        );
        
        % durations = seconds(intvls.EndTime - intvls.StartTime);
        % current_UE.SimData.TimeStats = struct(...
        %     'Mean',   mean(durations), ...
        %     'Max',    max(durations), ...
        %     'Min',    min(durations), ...
        %     'Total',  sum(durations), ...
        %     'Num_Passes', numel(durations) ... 
        % );
        
        UEs{idx} = current_UE;
        send(dq, []);
    end
    fprintf('\nGeometry calculation complete (%.1f sec).\n', toc);
    
    % Coverage Stats
    nT = length(UEs{1}.SimData.Time);
    counts = zeros(NumUEs, nT);

    for idx = 1:NumUEs
        counts(idx,:) = UEs{idx}.SimData.Num_visible;
    end

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
    metrics.Num_visible = counts; % returns num visible timeseries
    
    % Grabs 'SimData' from each cell and stores it in a new cell array
    metrics.SimData = cellfun(@(x) x.SimData, UEs, 'UniformOutput', false);
    % metrics.SimData = UEs.SimData;

    metrics.throughput_10pct  = 69;
    metrics.throughput_mean   = 420;


    %% Link Budget Calculation
    if calc_link
        tic
        parfor (idx = 1:NumUEs,num_workers)
            UEs{idx}.DL = link_calc(UEs{idx}, Cfg.DL, Cfg);
            send(dq, []);
        end
        fprintf('\nLoss calculation complete (%.1f sec).\n', toc);
    
        %% Data Aggregation & Metrics
        all_snr  = [];
        all_thpt = [];
        all_loss = [];
        all_adjusted_power_dBm = [];
        all_el_deg = [];
        all_pfd_W_MHz = [];
        
        
    
        meanThroughput = zeros(1,NumUEs);
        
        for idx = 1:NumUEs
            all_snr  = [all_snr, [UEs{idx}.DL.SNR]];
            all_thpt = [all_thpt, [UEs{idx}.DL.Throughput]];
            meanThroughput(idx) = mean(UEs{idx}.DL.Throughput);
            all_loss = [all_loss, [UEs{idx}.DL.Total_loss]];
            all_adjusted_power_dBm = [all_adjusted_power_dBm, [UEs{idx}.DL.Adjusted_EIRP_dBm]];
            all_el_deg =  [all_el_deg, [UEs{idx}.SimData.Elevation_deg]];
            all_pfd_W_MHz = [all_pfd_W_MHz, [UEs{idx}.DL.PFD_W_MHz]];
        end
    

        metrics.throughput_10pct  = prctile(all_thpt, 10);
        metrics.throughput_mean   = mean(all_thpt);
    end

    %% Plot Generation & Saving
    if plot_results
        if (isfield(Cfg, 'Save_dir')) %Provide location to save results
            out_dir = Cfg.Save_dir;
        else
            %% Create Output Directory
            % Format: simulation_output/orbit_sats_inclination_phasing_date
            date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            folder_name = sprintf('%.0f_%d_%.0f_%.1f_%s', ...
                Cfg.Orbit_height/1e3, Cfg.Total_sats, Cfg.Inclination, Cfg.Phasing, date_str);
            out_dir = fullfile('simulation_output', folder_name);
            
            if ~exist(out_dir, 'dir')
                mkdir(out_dir);
            end
        end
    
        % Save the configuration file immediately
        save(fullfile(out_dir, 'Config.mat'), 'Cfg');



        num_plots = 10; %plot progress bar
        plot_dq = parallel.pool.DataQueue;
        updateLiveScriptProgress(num_plots, true); 
        afterEach(plot_dq, @(~) updateLiveScriptProgress(num_plots, false));
        tic

        if Cfg.DL.B < 1e6
            thpt_scale = 1e3;
            thpt_unit = 'kbps';
            b_unit = 'kHz';
        else
            thpt_scale = 1e6;
            thpt_unit = 'Mbps';
            b_unit = 'MHz';
        end
        
        % Combined Constellation Stats
        f1 = figure('Visible', 'off', 'Name', 'Combined Constellation Stats', 'Color', 'w', 'Position', [100 100 1000 600]); 
        t2 = tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Tile 1: Throughput PDF
        nexttile; histogram(all_thpt./thpt_scale, 'Normalization', 'pdf', 'FaceColor', '#D95319', 'EdgeColor', 'none');
        grid on; title('Global Throughput PDF (DL)'); xlabel(sprintf('Throughput (%s)', thpt_unit)); ylabel('Probability Density');
        
        % Tile 2: Throughput CDF
        nexttile; [f_thpt, x_thpt] = ecdf(all_thpt./thpt_scale); plot(x_thpt, f_thpt, 'LineWidth', 2, 'Color', '#7E2F8E');
        grid on; title('Throughput CDF (DL)'); xlabel(sprintf('Throughput (%s)', thpt_unit)); ylabel('Probability \leq x'); xlim([0 max(x_thpt)]);
        
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
            ['Type:            ' walker_str];
            ['Total Sats:      ' num2str(Cfg.Total_sats)];
            ['Planes/Sats:     ' num2str(Cfg.Num_planes) ' / ' num2str(Cfg.Sats_per_plane)];
            ['Inclination:     ' num2str(Cfg.Inclination, '%.1f') '\circ'];
            ['Phasing:         ' num2str(Cfg.Phasing)];
            ['Altitude:        ' num2str(Cfg.Orbit_height/1e3, '%.0f') ' km'];
            ['Min Elev:        ' num2str(Cfg.Min_elevation_UE, '%.1f') '\circ'];
        };
        
        % Column 2: Link Budget / RF Data
        col2_str = {
            '\bfLink Budget Specs\rm';
            ['Direction:      ' char(Cfg.DL.Direction)];
            ['Freq / BW:      ' num2str(Cfg.DL.f/1e9, '%.2f') ' GHz / ' num2str(Cfg.DL.B/thpt_scale, '%.1f') sprintf(' %s', b_unit)];
            ['Tx Type/Gain:   ' char(Cfg.DL.Tx_type) '  / ' num2str(Cfg.DL.G_tx, '%.1f') ' dBi'];
            ['P\_tx / EIRP:    ' num2str(Cfg.DL.Max_P_tx_dBm, '%.1f') ' dBm / ' num2str(Cfg.DL.Max_EIRP_dBm, '%.1f') ' dBm'];
            ['Rx Type/Gain:   ' char(Cfg.DL.Rx_type) '  / ' num2str(Cfg.DL.G_rx, '%.1f') ' dBi'];
            ['Noise Fig:      ' num2str(Cfg.DL.NF, '%.1f') ' dB'];
        };  
        
        % Column 3: Performance Results
        col3_str = {
            '\bfPerformance Results\rm';
            ['  10% Throughput:   ' num2str(metrics.throughput_10pct / thpt_scale, '%.2f') sprintf(' %s', thpt_unit)];
            ['  Mean Throughput:  ' num2str(metrics.throughput_mean / thpt_scale, '%.2f') sprintf(' %s', thpt_unit)];
            ['  Worst Coverage:   ' num2str(metrics.worst_coverage_percent, '%.2f') ' %'];
            ['\bfSNR (dB)\rm'];
            ['  Mean: ' num2str(mean(all_snr,'omitnan'), '%.2f')];
            ['  Min:  ' num2str(min(all_snr), '%.2f')];
            ['  Max:  ' num2str(max(all_snr), '%.2f')];
        };
        
        % Plot text in 3 evenly spaced columns
        text(0.05, 0.95, col1_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        text(0.38, 0.95, col2_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        text(0.72, 0.95, col3_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
        
        exportgraphics(f1, fullfile(out_dir, 'Global_Stats.png'), 'Resolution', 300);
        close(f1);
        send(plot_dq, []);
    
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
            set(gca, 'FontSize', 14);
            exportgraphics(f2, fullfile(out_dir, 'Map_Min_Sats.png'), 'Resolution', 300);
            close(f2);
            send(plot_dq, []);
            
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
            set(gca, 'FontSize', 14);
            exportgraphics(f3, fullfile(out_dir, 'Map_Mean_Sats.png'), 'Resolution', 300);
            close(f3);
            send(plot_dq, []);
            
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
            set(gca, 'FontSize', 14);
            exportgraphics(f4, fullfile(out_dir, 'Map_Coverage_Prob.png'), 'Resolution', 300);
            close(f4);
            send(plot_dq, []);


            % %% Figure 5: Adjusted Power vs Elevation
            % f5 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 600]);
            % 
            % % 15 is the marker size. 'filled' and FaceAlpha=0.05 makes high-density areas pop out
            % scatter(all_adjusted_power_dBm - Cfg.DL.G_tx, all_el_deg, 15, 'filled', ...
            %     'MarkerFaceColor', '#0072BD', 'MarkerFaceAlpha', 0.1);
            % 
            % grid on; box on;
            % title('Elevation Angle vs. Adjusted Power', 'FontSize', 22, 'FontWeight', 'bold');
            % xlabel('Adjusted Power (dBm)', 'FontSize', 18);
            % ylabel('Elevation Angle (deg)', 'FontSize', 18);
            % 
            % % Thicken the axes and set font size
            % set(gca, 'FontSize', 14, 'LineWidth', 1.5);
            % exportgraphics(f5, fullfile(out_dir, 'adjusted_power.png'), 'Resolution', 300);
            % close(f5);
            % send(plot_dq, []);
            
            %% Figure 6: Adjusted Power vs SNR
            f6 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 600]);
            
            % Using an orange accent color for contrast
            scatter(all_adjusted_power_dBm-Cfg.DL.G_tx, all_snr, 15, 'filled', ...
                'MarkerFaceColor', '#D95319', 'MarkerFaceAlpha', 0.1);
            
            grid on; box on;
            title('SNR vs. Adjusted Power', 'FontSize', 22, 'FontWeight', 'bold');
            xlabel('Adjusted Power (dBm)', 'FontSize', 18);
            ylabel('SNR (dB)', 'FontSize', 18);
            
            set(gca, 'FontSize', 14, 'LineWidth', 1.5);
            exportgraphics(f6, fullfile(out_dir, 'adjusted_power_but_snr.png'), 'Resolution', 300);
            close(f6);
            send(plot_dq, []);

            %% Figure 7: PFD regulation
            f7 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 600]);
            
            scatter(all_el_deg, all_pfd_W_MHz, 15, 'filled', ...
                'MarkerFaceColor', '#D95319', 'MarkerFaceAlpha', 0.1);
            
            grid on; box on;
            title('Elevation vs. PFD', 'FontSize', 22, 'FontWeight', 'bold');
            xlabel('Elevation angle (deg)', 'FontSize', 18);
            ylabel('PFD (dBW/m^2/MHz)', 'FontSize', 18);
            ylim([all_pfd_W_MHz(1)-5, all_pfd_W_MHz(1)+5]);
            
            set(gca, 'FontSize', 14, 'LineWidth', 1.5);
            exportgraphics(f7, fullfile(out_dir, 'pfd_regulation.png'), 'Resolution', 300);
            close(f7);
            send(plot_dq, []);

            %% MAP 8: Mean Throughput
            % 1. Force everything to be a column vector (using :) to prevent mismatches
            lats = lat_vector(:);
            lons = lon_vector(:);
            thpt_vals = meanThroughput(:) / thpt_scale;
            
            % 2. Run griddata with consistent dimensions
            meanThroughputGrid = griddata(lats, lons, thpt_vals, LatG, LonG, 'cubic'); 
            
            f8 = figure('Visible', 'off', 'Color', 'w');
            ax3 = axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, ...
                        'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
            axis off;  
            
            geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
            surfm(LatG, LonG, meanThroughputGrid, 'FaceAlpha', 0.5);
            
            for i = 1:length(lats)
                textm(lats(i), lons(i), sprintf('%.1f', thpt_vals(i)), ...
                      'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                      'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
            end
            
            cb3 = colorbar; 
            % 3. Use a comma for caxis and [] for string concatenation in ylabel
            caxis([0, max([thpt_vals; 1])]); 
            ylabel(cb3, ['Mean Throughput (', thpt_unit, ')']); 
            
            set(gca, 'FontSize', 14);
            exportgraphics(f8, fullfile(out_dir, 'Map_Mean_Throughput.png'), 'Resolution', 300);
            close(f8);
            send(plot_dq, []);

            
            %% MAP 9: Elevation distribution
            f9 = figure('Visible', 'off', 'Color', 'w');
            histogram(all_el_deg);
            grid on; box on;
            xlabel('Elevation Angle (deg)', 'FontWeight', 'bold');
            ylabel('Number of Occurences', 'FontWeight', 'bold');
            title('Distribution of Elevation Angles', 'FontSize', 14);

            exportgraphics(f9, fullfile(out_dir, 'Elevation_distribution.png'), 'Resolution', 300);
            close(f9);
            send(plot_dq, []);

            %% MAP 10: Elevation distribution
            % Calculate nadir steering angles
            Re = 6378.14e3;
            r = Re + Cfg.Orbit_height; % Orbit radius
            % Nadir Angle (eta) - The tilt of the satellite antenna
            eta = asin((Re/r) * cos(all_el_deg));

            f10 = figure('Visible', 'off', 'Color', 'w');
            histogram(eta);
            grid on; box on;
            xlabel('Nadir Steering Angle (deg)', 'FontWeight', 'bold');
            ylabel('Number of Occurences', 'FontWeight', 'bold');
            title('Distribution of Nadir Angles', 'FontSize', 14);

            exportgraphics(f10, fullfile(out_dir, 'nadir_steering_distribution.png'), 'Resolution', 300);
            close(f10);
            send(plot_dq, []);


            
        catch ME
            fprintf('Warning: Map generation failed. Error message:\n');
            disp(ME.message);
        end
    
        fprintf('\nPlots generation complete (%.1f sec).', toc);
    end
end


