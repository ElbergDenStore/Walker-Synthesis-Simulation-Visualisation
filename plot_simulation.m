function plot_simulation(metrics, use_parallel)
    Cfg = metrics.Cfg;
    UEs = metrics.UEs;

    %% 2. Setup Output Directory
    if isfield(Cfg, 'Save_dir')
        out_dir = Cfg.Save_dir;
    else
        date_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        folder_name = sprintf('Plots_%.0fkm_%d_Sats_%s', Cfg.Orbit_height/1e3, Cfg.Total_sats, date_str);
        out_dir = fullfile('simulation_output', folder_name);
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    end
    save(fullfile(out_dir, 'Config.mat'), 'Cfg');

    %% 3. Extract Arrays from UEs Struct (Fast Vectorization)
    fprintf('Extracting plotting data from metrics...\n');
    lat_vector = [UEs.Lat]';
    lon_vector = [UEs.Lon]';
    
    SimDataArray = [UEs.SimData];
    all_el_deg   = single(vertcat(SimDataArray.Elevation_deg));
    
    % Determine Throughput Scaling
    if Cfg.DL.B < 1e6
        thpt_scale = 1e3; thpt_unit = 'kbps'; b_unit = 'kHz';
    else
        thpt_scale = 1e6; thpt_unit = 'Mbps'; b_unit = 'MHz';
    end

    % Extract Link Data if available
    if isfield(UEs(1), 'DL')
        DL_Array = [UEs.DL];
        all_thpt = vertcat(DL_Array.Throughput);
        all_snr  = vertcat(DL_Array.SNR);
        all_sinr = vertcat(DL_Array.SINR);
        
        valid_thpt = all_thpt(~isnan(all_thpt));
        meanThroughput = mean(all_thpt, 2, 'omitnan');
        
        % Calculate 10% and Mean for stats
        metrics.throughput_10pct = prctile(valid_thpt, 10);
        metrics.throughput_mean  = mean(valid_thpt);
    end

    %% 4. Pre-calculate Shared Map Data
    lat_lim = [min(Cfg.Flat_UE_array.Lats) max(Cfg.Flat_UE_array.Lats)];
    lon_lim = [min(Cfg.Flat_UE_array.Lons) max(Cfg.Flat_UE_array.Lons)];
    nLat = 500; nLon = 500; 
    [LonG, LatG] = meshgrid(linspace(lon_lim(1), lon_lim(2), nLon), linspace(lat_lim(1), lat_lim(2), nLat));
    
    % Read shapefile once and broadcast to workers to prevent file lock errors
    land = shaperead('landareas.shp', 'UseGeoCoords', true);

    %% 5. Parallel Plot Generation
    num_tasks = 7;
    fprintf('Generating %d plots in parallel...\n', num_tasks);
    
    % Set up progress bar for parallel pool
    updateLiveScriptProgress(num_tasks, true);
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(~) updateLiveScriptProgress(num_tasks, false));
    
    if use_parallel
        num_workers = Inf; 
    else
        num_workers = 0;   
    end
    


    tic;
    
    parfor(task_id = 1:num_tasks, num_workers)
        try
            switch task_id
                case 1
                    generate_global_stats(Cfg, metrics, all_thpt, all_snr, all_sinr, thpt_scale, thpt_unit, b_unit, out_dir);
                case 2
                    generate_map_min_sats(metrics.minNumberSatellites, lat_vector, lon_vector, LonG, LatG, lat_lim, lon_lim, land, Cfg.NumUEs, out_dir);
                case 3
                    generate_map_mean_sats(metrics.meanNumberSatellites, lat_vector, lon_vector, LonG, LatG, lat_lim, lon_lim, land, Cfg.NumUEs, out_dir);
                case 4
                    generate_map_coverage(metrics.prob_coverage, lat_vector, lon_vector, LonG, LatG, lat_lim, lon_lim, land, Cfg.NumUEs, out_dir);
                case 5
                    generate_map_throughput(meanThroughput, lat_vector, lon_vector, LonG, LatG, lat_lim, lon_lim, land, Cfg.NumUEs, thpt_scale, thpt_unit, out_dir);
                case 6
                    generate_elevation_dist(all_el_deg, out_dir);
                case 7
                    generate_nadir_dist(all_el_deg, Cfg.Orbit_height, out_dir);
            end
        catch ME
            fprintf('Task %d failed: %s\n', task_id, ME.message);
        end
        send(dq, []); % Update progress bar
    end
    fprintf('Plots generation complete (%.1f sec).\n', toc);
end


%% ========================================================================
%% HELPER FUNCTIONS (Local Functions)
%% ========================================================================

function generate_global_stats(Cfg, metrics, all_thpt, all_snr, all_sinr, thpt_scale, thpt_unit, b_unit, out_dir)
    f1 = figure('Visible', 'off', 'Name', 'Combined Constellation Stats', 'Color', 'w', 'Position', [100 100 1000 600]); 
    t2 = tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    nexttile; histogram(all_thpt./thpt_scale, 'Normalization', 'pdf', 'FaceColor', '#D95319', 'EdgeColor', 'none');
    grid on; title('Global Throughput PDF (DL)'); xlabel(sprintf('Throughput (%s)', thpt_unit)); ylabel('Probability Density');
    
    nexttile; [f_thpt, x_thpt] = ecdf(all_thpt(:)./thpt_scale); plot(x_thpt, f_thpt, 'LineWidth', 2, 'Color', '#7E2F8E');
    grid on; title('Throughput CDF (DL)'); xlabel(sprintf('Throughput (%s)', thpt_unit)); ylabel('Probability \leq x'); 
    
    nexttile(3, [1 2]); axis off; 
    
    walker_str = 'Walker Delta';
    if Cfg.WalkerStar; walker_str = 'Walker Star'; end

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
    
    col2_str = {
        '\bfLink Budget Specs\rm';
        ['Direction:      ' char(Cfg.DL.Direction)];
        ['Freq / BW:      ' num2str(Cfg.DL.f/1e9, '%.2f') ' GHz / ' num2str(Cfg.DL.B/thpt_scale, '%.1f') sprintf(' %s', b_unit)];
        ['Tx Type/Gain:   ' char(Cfg.DL.Tx_type) '  / ' num2str(Cfg.DL.G_tx, '%.1f') ' dBi'];
        ['P\_tx / EIRP:    ' num2str(Cfg.DL.Max_P_tx_dBm, '%.1f') ' dBm / ' num2str(Cfg.DL.Max_EIRP_dBm, '%.1f') ' dBm'];
        ['Rx Type/Gain:   ' char(Cfg.DL.Rx_type) '  / ' num2str(Cfg.DL.G_rx, '%.1f') ' dBi'];
        ['Noise Fig:      ' num2str(Cfg.DL.NF, '%.1f') ' dB'];
        ['RU / FRF:       ' num2str(Cfg.RU*100, '%.0f') '% / ' num2str(Cfg.FRF)];
    };  
    
    valid_snr  = all_snr(~isnan(all_snr));
    valid_sinr = all_sinr(~isnan(all_sinr));
    
    col3_str = {
        '\bfPerformance Results\rm';
        ['  10% Throughput:   ' num2str(metrics.throughput_10pct / thpt_scale, '%.2f') sprintf(' %s', thpt_unit)];
        ['  Mean Throughput:  ' num2str(metrics.throughput_mean / thpt_scale, '%.2f') sprintf(' %s', thpt_unit)];
        ['  Worst Coverage:   ' num2str(metrics.worst_coverage_percent, '%.2f') ' %'];
        ['\bfSNR / SINR (dB)\rm'];
        ['  Mean: ' num2str(mean(valid_snr), '%.2f') '  /  ' num2str(mean(valid_sinr), '%.2f')];
        ['  Min:  ' num2str(min(valid_snr), '%.2f') '  /  ' num2str(min(valid_sinr), '%.2f')];
        ['  Max:  ' num2str(max(valid_snr), '%.2f') '  /  ' num2str(max(valid_sinr), '%.2f')];
    };
    
    text(0.05, 0.95, col1_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
    text(0.38, 0.95, col2_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
    text(0.72, 0.95, col3_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Consolas', 'Interpreter', 'tex');
    
    exportgraphics(f1, fullfile(out_dir, 'Global_Stats.png'), 'Resolution', 300);
    close(f1);
end

function generate_map_min_sats(vals, lat_v, lon_v, LonG, LatG, lat_lim, lon_lim, land, nUEs, out_dir)
    ValG = griddata(lon_v, lat_v, vals, LonG, LatG, 'cubic');
    f = figure('Visible', 'off', 'Color', 'w');
    axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, 'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
    axis off;  
    geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
    surfm(LatG, LonG, ValG, 'FaceAlpha', 0.5);
    if nUEs < 50
        for i = 1:nUEs
            textm(lat_v(i), lon_v(i), sprintf('%d', vals(i)), 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k'); 
        end
    end
    cb = colorbar; caxis([0 max([vals(:); 1])]); ylabel(cb, 'Min Number of Satellites');
    set(gca, 'FontSize', 14);
    exportgraphics(f, fullfile(out_dir, 'Map_Min_Sats.png'), 'Resolution', 300);
    close(f);
end

function generate_map_mean_sats(vals, lat_v, lon_v, LonG, LatG, lat_lim, lon_lim, land, nUEs, out_dir)
    ValG = griddata(lon_v, lat_v, vals, LonG, LatG, 'cubic');
    f = figure('Visible', 'off', 'Color', 'w');
    axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, 'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
    axis off;  
    geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
    surfm(LatG, LonG, ValG, 'FaceAlpha', 0.5);
    if nUEs < 50
        for i = 1:nUEs
            textm(lat_v(i), lon_v(i), sprintf('%.1f', vals(i)), 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
    cb = colorbar; caxis([0 max([vals(:); 1])]); ylabel(cb, 'Mean Number of Satellites');
    set(gca, 'FontSize', 14);
    exportgraphics(f, fullfile(out_dir, 'Map_Mean_Sats.png'), 'Resolution', 300);
    close(f);
end

function generate_map_coverage(vals, lat_v, lon_v, LonG, LatG, lat_lim, lon_lim, land, nUEs, out_dir)
    ValG = griddata(lon_v, lat_v, vals, LonG, LatG, 'cubic');
    f = figure('Visible', 'off', 'Color', 'w');
    axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, 'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
    axis off; 
    geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
    surfm(LatG, LonG, ValG, 'FaceAlpha', 0.5);
    if nUEs < 50
        for i = 1:nUEs
            textm(lat_v(i), lon_v(i), sprintf('%.0f', vals(i)), 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k'); 
        end
    end
    cb = colorbar; caxis([0 100]); ylabel(cb, 'P(Sat \geq 1) [%]');
    set(gca, 'FontSize', 14);
    exportgraphics(f, fullfile(out_dir, 'Map_Coverage_Prob.png'), 'Resolution', 300);
    close(f);
end

function generate_map_throughput(vals, lat_v, lon_v, LonG, LatG, lat_lim, lon_lim, land, nUEs, thpt_scale, thpt_unit, out_dir)
    lats = lat_v(:); lons = lon_v(:); thpt_vals = vals(:) / thpt_scale;
    ValG = griddata(lats, lons, thpt_vals, LatG, LonG, 'cubic'); 
    
    [~, dist_to_nearest_UE] = dsearchn([lons, lats], [LonG(:), LatG(:)]);
    mask = reshape(dist_to_nearest_UE, size(LonG)) > 2.5;
    ValG(mask) = NaN; % choose between unlimited interpolation or masked
    
    f = figure('Visible', 'off', 'Color', 'w');
    axesm('lambertstd', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim, 'Frame', 'on', 'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on');
    axis off;  
    geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8]);
    surfm(LatG, LonG, ValG, 'FaceAlpha', 0.5);
    if nUEs < 50
        for i = 1:length(lats)
            textm(lats(i), lons(i), sprintf('%.1f', thpt_vals(i)), 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
    cb = colorbar; caxis([0, max([thpt_vals; 1])]); ylabel(cb, ['Mean Throughput (', thpt_unit, ')']); 
    set(gca, 'FontSize', 14);
    exportgraphics(f, fullfile(out_dir, 'Map_Mean_Throughput.png'), 'Resolution', 300);
    close(f);
end

function generate_elevation_dist(all_el_deg, out_dir)
    f = figure('Visible', 'off', 'Color', 'w');
    histogram(all_el_deg);
    grid on; box on;
    xlabel('Elevation Angle (deg)', 'FontWeight', 'bold');
    ylabel('Number of Occurrences', 'FontWeight', 'bold');
    title('Distribution of Elevation Angles', 'FontSize', 14);
    exportgraphics(f, fullfile(out_dir, 'Elevation_distribution.png'), 'Resolution', 300);
    close(f);
end

function generate_nadir_dist(all_el_deg, orbit_height, out_dir)
    Re = 6378.14e3;
    r = Re + orbit_height; 
    eta = rad2deg(asin((Re/r) * cosd(all_el_deg)));

    f = figure('Visible', 'off', 'Color', 'w');
    histogram(eta);
    grid on; box on;
    xlabel('Nadir Steering Angle (deg)', 'FontWeight', 'bold');
    ylabel('Number of Occurrences', 'FontWeight', 'bold');
    title('Distribution of Nadir Angles', 'FontSize', 14);
    exportgraphics(f, fullfile(out_dir, 'nadir_steering_distribution.png'), 'Resolution', 300);
    close(f);
end