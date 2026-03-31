% %% 1. Configuration & Region Selection
% % Toggle between 'Nordjylland', 'Denmark', or 'Full'
REGION = 'Full3000'; 

switch REGION
    case 'Nordjylland'
        latlim = [56.5, 58.0];
        lonlim = [8.0, 11.0];
        people_per_ue = 300;
        dataFile = 'Nordjylland300.mat';
    case 'Denmark'
        latlim = [54.5, 58.0];
        lonlim = [8.0, 15.5];
        people_per_ue = 3000;
        dataFile = 'Denmark3000.mat';
    case 'Full300'
        latlim = [54.5, 83.9]; 
        lonlim = [-60, 30.0];
        people_per_ue = 300;
        dataFile = 'Full300.mat';
    case 'Full3000'
        latlim = [54.5, 83.9]; 
        lonlim = [-60, 30.0];
        people_per_ue = 3000;
        dataFile = 'Full3000.mat';
    case 'Full30000'
        latlim = [54.5, 83.9]; 
        lonlim = [-60, 30.0];
        people_per_ue = 30000;
        dataFile = 'Full30000.mat';
end

FORCE_RERUN = true; 

if exist(dataFile, 'file') && ~FORCE_RERUN 
    fprintf('Loading %s data from %s...\n', REGION, dataFile);
    load(dataFile);
else
    fprintf('Running new simulation for %s...\n', REGION);
    constellation_type = "walkerDelta";
    ue_grid_size = "medium";
    duration = "medium";
    frequency = "ku";
    height_km = 1000;
    Cfg = get_cfg(height_km,constellation_type,ue_grid_size,duration,frequency);

    
    Cfg.Accept_Flat_UE_array = true; Cfg.Equal_UE_area = false; 
    [Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, Total_Pop] = generate_population_based_UEs(latlim, lonlim, people_per_ue);


    % VIP_UEs = aalborg, nuuk, kastrup? anden mili base? a fleet of ships
    % Cfg.Flat_UE_array.Lats = [Cfg.Flat_UE_array.Lats; VIP_UEs.Lats]
    % Cfg.Flat_UE_array.Lons = [Cfg.Flat_UE_array.Lons; VIP_UEs.Lons]
    calc_link = true; 
    Cfg.Use_P618 = false; %simple atmospheric loss
    Cfg.Simple_Atmospheric_Loss_dB = 1;
    Cfg.FRF = 3;
    Cfg.RU = 1;
    Cfg.Share_bandwidth  = true; 
    plot_results = true; 
    use_parallel = true;
    metrics = coverage_simulator_function(Cfg, use_parallel, calc_link);

    plot_simulation(metrics, use_parallel);
    % save(dataFile);
end

%%%%% VISUALISE UES
% --- VISUAL VERIFICATION PLOT ---
fprintf('Plotting UE distribution map...\n');

% Create a clean, white-background figure window
figure('Name',  'UE Distribution', 'Color', 'w', 'Position', [100, 100, 800, 800]);

% Create geographic axes
gx = geoaxes;

% Scatter the UEs: Size 30, Red, Filled, with a slight black edge for visibility
geoscatter(gx, Cfg.Flat_UE_array.Lats, Cfg.Flat_UE_array.Lons, 30, 'r', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Lock the map view specifically to Greenland's coordinates
geolimits(gx, latlim, lonlim);

% Add a basemap to show the terrain and ice sheet 
% (Other good options: 'topographic', 'streets-light', or 'bluegreen')
geobasemap(gx, 'satellite'); 

% Add a dynamic title
Cfg.NumUEs = length(Cfg.Flat_UE_array.Lats);
title(sprintf(' UE Distribution\nTotal UEs: %d (1 per %d people)', ...
    Cfg.NumUEs, round(Total_Pop/Cfg.NumUEs)), 'FontWeight', 'bold', 'FontSize', 14);

pause(3*(Cfg.NumUEs/3000)) %to load UE distribution, 3000UEs take 3 seconds
exportgraphics(gcf, 'screenshots/UEDistribution.png', 'Resolution', 300);