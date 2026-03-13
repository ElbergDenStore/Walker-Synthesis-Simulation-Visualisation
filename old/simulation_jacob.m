%% Constellation Design Tool

clear all;
close all;
clc;

fprintf('Defining constellation...\n')
sc = satelliteScenario;

sc.StartTime = datetime('1-Jun-2025 00:00:00', 'TimeZone', 'UTC');
sc.StopTime  = datetime('1-Jun-2025 23:59:59', 'TimeZone', 'UTC');

% Sample Time --> Make it smaller for more accuracy (also slower execution time)
sc.SampleTime = 60;
step = sc.SampleTime;

%Elevation angle --> It adjust the footprint area. Somewhere between 20 and
%40 is reasonable. The lower it is, the larger the footprint
minElevationUE = 20;

%Example of Constellation Definition
orbit_height = 700e3;
r_earth = 6378.14e3;
% Walker-Delta Constellation: 64 satellites, 8 planos, phasing = 2, i = 75°
%sat = walkerStar(scenario,radius            ,  inclination, totalSatellites   , geometryPlanes, phasing)
%sat = walkerDelta(sc        ,1000e3+r_earth  , 75         , 4                , 4             , 0       ,Name="Sat", OrbitPropagator="sgp4");
sat = walkerDelta(sc       ,orbit_height+r_earth  , 75         , 64               , 8             , 3       ,Name="S4D", OrbitPropagator="sgp4"); 

%Example of Code for WalkerStar Definition
%sat = walkerStar(sc      ,1000e3+6378.14e3  , 70         , 64               , 8             , 0       ,Name="S4D", OrbitPropagator="sgp4");


%% Constellation Visualization
%WARNING: If you run the next line, make sure that the code only run from
%line 1-36. Otherwise, MatLab tend to crash when adding UE/GS to the
%visualization tool (few lines below)
%satelliteScenarioViewer(sc);

%% UE/GS Grid

%With linspace you choose the range of latitudes and longitudes and how
%many points you want to visualize in the grid. They can emulate User
%Equipment or Ground Stations

%Right now, latitudes and longitudes corresponds to Denmark
lat_vec = linspace(50, 85, 1);       % 6 latitudes
lon_vec = linspace(-60, 30, 1);       % 6 longitudes
[LonGrid, LatGrid] = meshgrid(lon_vec, lat_vec);
nPts = numel(LonGrid);

%From this point, you should not change the code unless you want to
%generate new results. Otherwise, simulation parameters are above this
%point.

%% Coverage Simulation

fprintf('Defining UEs and links to the satellite...\n')

SatData = cell(1, nPts); %preallocating, might be faster?

for pt = 1:nPts
    lat = LatGrid(pt)
    lon = LonGrid(pt)
    ue = groundStation(sc,lat,lon, 'Name',['UE_', num2str(pt)], 'MinElevationAngle', minElevationUE)
    % We calculate the links between the satelites and user
    % equipment/ground station
    ac = access(sat, ue)
    %Save intervals in which UE-SAT links establish connectivity
    intvls = accessIntervals(ac)
    % intvl{pt} = intvls_temp
    

    %%%%% ADDED LOGIC FOR AER (azimuth, elevation and range)
    ue_results = struct('Time', {}, 'SatID', {}, 'Range_km', {}, 'Elevation_deg', {});
    count = 0;
    
    if isempty(intvls)
         SatData{pt} = [];
    else
        % 1. EXACT CALCULATION
        raw_steps = seconds(intvls.EndTime - intvls.StartTime) ./ sc.SampleTime;
        num_steps = round(raw_steps); 
        n_exact   = sum(num_steps + 1);
        
        % 2. PRE-ALLOCATE VECTORS
        vec_Time      = NaT(1, n_exact); 
        vec_SatID     = zeros(1, n_exact);
        vec_Range     = zeros(1, n_exact);
        vec_Elevation = zeros(1, n_exact);
        
        count = 0;
        
        for row = 1:height(intvls)
            % Extract Sat ID
            sourceName = string(intvls.Source(row)); 
            sat_idx = sscanf(sourceName, "S4D_%d");
            
            t_start = intvls.StartTime(row);
            t_end   = intvls.EndTime(row);
            
            t_idx = find(times >= t_start & times <= t_end);
            
            if ~isempty(t_idx)
                relevant_times = times(t_idx);
                num_new_points = numel(relevant_times);
                
                for t_k = 1:num_new_points
                    t_scalar = relevant_times(t_k);
                    
                    [~, el, r] = aer(ue, sat(sat_idx), t_scalar);
                    
                    count = count + 1;
                    
                    vec_Time(count)      = t_scalar;
                    vec_SatID(count)     = sat_idx;
                    vec_Range(count)     = r / 1000;
                    vec_Elevation(count) = el;
                end
            end
        end
        
        % 3. STORE AS SINGLE STRUCT OF ARRAYS (Fastest)
        if count > 0
            % We truncate the vectors to 1:count in case of any estimation mismatch
            SatData{pt} = struct(...
                'Time',          vec_Time(1:count), ...
                'SatID',         vec_SatID(1:count), ...
                'Range_km',      vec_Range(1:count), ...
                'Elevation_deg', vec_Elevation(1:count) ...
            );
        else
            SatData{pt} = [];
        end
    end
    
    % Save to main list
    SatData{pt} = ue_results;
    fprintf('Progress: %.2f%%...\n', 100*pt/nPts);


end

% AI generated statistics
%% ---------------------------------------------------------
%  STATISTICS & VISUALIZATION MODULE
%  ---------------------------------------------------------
fprintf('Generating Statistics...\n');

% 1. Data Extraction
% We need to pull the numbers out of the 'SatData' structures
all_ranges_km = [];
all_elevs_deg = [];
ue_legend_names = {};

% Create a figure for the "Per UE" analysis
figure('Name', 'Per-UE Distributions', 'Color', 'w');
tiledlayout(2,1);

% --- PLOT 1: SLANT RANGE (Per UE) ---
nexttile;
hold on; grid on;
title('Slant Range Distribution (Per UE)');
xlabel('Slant Range (km)'); 
ylabel('Frequency (Counts)');

% Iterate through UEs to plot individual histograms
for pt = 1:nPts
    % Check if this UE has data
    if ~isempty(SatData{pt})
        % Extract vector of ranges from the struct array
        % [struct.Field] creates a standard vector
        r_vec = [SatData{pt}.Range_km];
        
        % Plot histogram (Adjust 'BinWidth' as needed)
        % 'FaceAlpha' makes it transparent so you can see overlaps
        histogram(r_vec, 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', ['UE_', num2str(pt)]);
        
        % Collect for the "Combined" plot later
        all_ranges_km = [all_ranges_km, r_vec];
        ue_legend_names{end+1} = ['UE_', num2str(pt)];
    end
end
legend(ue_legend_names, 'Location', 'bestoutside');
hold off;

% --- PLOT 2: ELEVATION ANGLE (Per UE) ---
nexttile;
hold on; grid on;
title('Elevation Angle Distribution (Per UE)');
xlabel('Elevation Angle (deg)'); 
ylabel('Frequency (Counts)');

for pt = 1:nPts
    if ~isempty(SatData{pt})
        el_vec = [SatData{pt}.Elevation_deg];
        
        % Plot
        histogram(el_vec, 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', ['UE_', num2str(pt)]);
        
        % Collect for combined
        all_elevs_deg = [all_elevs_deg, el_vec];
    end
end
legend(ue_legend_names, 'Location', 'bestoutside');
hold off;


%% --- COMBINED (AGGREGATE) PLOTS ---
% This shows the "Global" performance (Averaged over all indices)

figure('Name', 'Combined Constellation Statistics', 'Color', 'w');
t = tiledlayout(2,2);

% 1. Combined Range Histogram
nexttile;
histogram(all_ranges_km, 'Normalization', 'pdf', 'FaceColor', '#0072BD', 'EdgeColor', 'none');
grid on;
title('Global Slant Range PDF');
xlabel('Range (km)'); ylabel('Probability Density');
xline(mean(all_ranges_km), '--r', ['Mean: ' num2str(round(mean(all_ranges_km))) 'km'], 'LineWidth', 1.5);

% 2. Combined Elevation Histogram
nexttile;
histogram(all_elevs_deg, 'Normalization', 'pdf', 'FaceColor', '#D95319', 'EdgeColor', 'none');
grid on;
title('Global Elevation PDF');
xlabel('Elevation (deg)'); ylabel('Probability Density');
xline(mean(all_elevs_deg), '--r', ['Mean: ' num2str(round(mean(all_elevs_deg))) 'deg'], 'LineWidth', 1.5);

% 3. Combined Basic Statistics Table (Printed on Plot)
nexttile([1,2]); % Span 2 columns
axis off;

% Calculate stats
stats_str = {
    ['\bfConstellation Statistics Summary\rm'];
    ['--------------------------------'];
    ['\bfTotal Access Events:\rm ' num2str(numel(all_ranges_km))];
    [''];
    ['\bfSlant Range (km):\rm'];
    ['  Mean: ' num2str(mean(all_ranges_km), '%.2f')];
    ['  Min:  ' num2str(min(all_ranges_km), '%.2f')];
    ['  Max:  ' num2str(max(all_ranges_km), '%.2f')];
    [''];
    ['\bfElevation (deg):\rm'];
    ['  Mean: ' num2str(mean(all_elevs_deg), '%.2f')];
    ['  Min:  ' num2str(min(all_elevs_deg), '%.2f')];
    ['  Max:  ' num2str(max(all_elevs_deg), '%.2f')];
};

text(0.5, 0.5, stats_str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Interpreter', 'tex');


 % % Simulation times
 % times = sc.StartTime:seconds(step):sc.StopTime;
 % 
 % %We iterate to see access times
 % nT = numel(times);
 % counts = zeros(nPts,nT);
 % 
 % for pt = 1: nPts
 %    intvls_temp = intvl{pt};
 % 
 %    for k = 1:nT
 %    t = times(k);
 %    active = (intvls_temp.StartTime <= t) & (intvls_temp.EndTime >= t);
 %    counts(pt,k) = sum(active);
 %    end
 % end


 %% Results

%Minimum number of satellites
minNumberSatellites = min(counts,[], 2);

%Average number of satellites
meanNumberSatellites = mean(counts, 2);

%Probability at least >= 1 satelite

more_1_satellites = counts >= 1;
more_1_satellites = 100*sum(more_1_satellites,2)./nT;

%Max Gap w/o coverage

maxGapMinutes = zeros(nPts,1);
for pt = 1:nPts
    row = counts(pt,:) == 0;  
    maxZeroStreak = 0;
    currentStreak = 0;
    
    for k = 1:length(row)
        if row(k)
            currentStreak = currentStreak + 1;
            maxZeroStreak = max(maxZeroStreak, currentStreak);
        else
            currentStreak = 0;
        end
    end
    
    maxGapMinutes(pt) = maxZeroStreak * step / 60;
end

%% Definition of TLE 

% Initialise:
N = length(sat);
TLE = cell(N,1);
% Loop over satellites to extract TLEs (for Matlab):
for j = 1 : N
    TLEmatlab{j} = getTLEmatlab(sc.Satellites(j));
end

%% Definition of vector for latitute and longitude

 for pt = 1:nPts
     lat_vector(pt) = LatGrid(pt);
     lon_vector(pt) = LonGrid(pt);
 end

%% Figure --> Min Numero Satelites (Surf)

%1) Geographical Grid
lat_lim = [min(lat_vec) max(lat_vec)];
lon_lim = [min(lon_vec) max(lon_vec)];
nLat = 200;  % Grid resolution
nLon = 200;  % Grid resolution

[LonG, LatG] = meshgrid( linspace(lon_lim(1), lon_lim(2), nLon) , linspace(lat_lim(1), lat_lim(2), nLat) );

%2) griddata interpolation
ValG = griddata(lon_vector, lat_vector, minNumberSatellites, LonG, LatG, 'cubic');

%3) Map draw
figure;
ax = axesm('lambertstd', ...
           'MapLatLimit', lat_lim, ...
           'MapLonLimit', lon_lim, ...
           'Frame', 'on', 'Grid', 'on', ...
           'MeridianLabel','on','ParallelLabel','on');
axis off;  

land = shaperead('landareas.shp','UseGeoCoords',true);
geoshow([land.Lat], [land.Lon], ...
        'DisplayType','polygon', ...
        'FaceColor',[0.8 0.8 0.8]);

% 4) heatmap
surfm(LatG, LonG, ValG,'FaceAlpha',0.5);

for pt = 1:nPts
    textm(lat_vector(pt), lon_vector(pt), ...
          sprintf('%d', minNumberSatellites(pt)), ...
          'HorizontalAlignment','center', ...
          'VerticalAlignment','middle', ...
          'FontSize', 14, ...
          'FontWeight', 'bold', ...
          'Color', 'k'); 
end

% 5) Colorbar and Title
cb = colorbar;
caxis([0 max([minNumberSatellites; 1])]);
ylabel(cb, 'Min Number of Satellites');
%title('Heatmap - Lambert Conformal Projection');
set(gca,'FontSize', 25);

%% Figure --> Mean Numero Satelites (Surf)

%1) Geographical Grid
lat_lim = [min(lat_vec) max(lat_vec)];
lon_lim = [min(lon_vec) max(lon_vec)];
nLat = 200;  % Grid Resolution
nLon = 200;  % Grid Resolution

[LonG, LatG] = meshgrid( linspace(lon_lim(1), lon_lim(2), nLon) , linspace(lat_lim(1), lat_lim(2), nLat) );

%2)  griddata interpolation
ValG = griddata(lon_vector, lat_vector, meanNumberSatellites, LonG, LatG, 'cubic');

%3) Map Draw
figure;
ax = axesm('lambertstd', ...
           'MapLatLimit', lat_lim, ...
           'MapLonLimit', lon_lim, ...
           'Frame', 'on', 'Grid', 'on', ...
           'MeridianLabel','on','ParallelLabel','on');
axis off;  

land = shaperead('landareas.shp','UseGeoCoords',true);
geoshow([land.Lat], [land.Lon], ...
        'DisplayType','polygon', ...
        'FaceColor',[0.8 0.8 0.8]);

% 4) heatmap
surfm(LatG, LonG, ValG,'FaceAlpha',0.5);

for pt = 1:nPts
    textm(lat_vector(pt), lon_vector(pt), ...
          sprintf('%.1f', meanNumberSatellites(pt)), ...
          'HorizontalAlignment','center', ...
          'VerticalAlignment','middle', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'Color', 'k');
end

% 5) Colorbar and title
cb = colorbar;
caxis([0 max(meanNumberSatellites)]);
ylabel(cb, 'Mean Number of Satellites');
%title('Heatmap - Lambert Conformal Projection');
set(gca,'FontSize', 25);

 %% Figure percentage >) 1 Satellite (Surf)
 
%1) Geographical Grid
lat_lim = [min(lat_vec) max(lat_vec)];
lon_lim = [min(lon_vec) max(lon_vec)];
nLat = 200;  % Grid Resolution
nLon = 200;  % Grid Resolution

[LonG, LatG] = meshgrid( linspace(lon_lim(1), lon_lim(2), nLon) , linspace(lat_lim(1), lat_lim(2), nLat) );

%2) griddata interpolation
ValG = griddata(lon_vector, lat_vector, more_1_satellites, LonG, LatG, 'cubic');

%3) Map Draw
figure;
ax = axesm('lambertstd', ...
            'FontSize', 20, ...
           'MapLatLimit', lat_lim, ...
           'MapLonLimit', lon_lim, ...
           'Frame', 'on', 'Grid', 'on', ...
           'MeridianLabel','on','ParallelLabel','on');
axis off; 

land = shaperead('landareas.shp','UseGeoCoords',true);
geoshow([land.Lat], [land.Lon], ...
        'DisplayType','polygon', ...
        'FaceColor',[0.8 0.8 0.8]);

% 4) Heatmap
surfm(LatG, LonG, ValG,'FaceAlpha',0.5);

for pt = 1:nPts
    textm(lat_vector(pt), lon_vector(pt), ...
          sprintf('%.f', more_1_satellites(pt)), ...
          'HorizontalAlignment','center', ...
          'VerticalAlignment','middle', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'Color', 'k'); 
end

% 5) Colorbar and Title
cb = colorbar;
caxis([0 100]);
ylabel(cb, 'P(Sat \geq 1)');
%title('Heatmap - Lambert Conformal Projection');
set(gca,'FontSize', 25);


%% Max Gap (Minutos) Sin Cobertura (Surf)

%1) Geographical Grid
lat_lim = [min(lat_vec) max(lat_vec)];
lon_lim = [min(lon_vec) max(lon_vec)];
nLat = 200;  % Grid Resolution
nLon = 200;  % Grid Resolution

[LonG, LatG] = meshgrid( linspace(lon_lim(1), lon_lim(2), nLon) , linspace(lat_lim(1), lat_lim(2), nLat) );

%2) griddata interpolation
ValG = griddata(lon_vector, lat_vector, maxGapMinutes, LonG, LatG, 'cubic');

%3) Map Draw
figure;
ax = axesm('lambertstd', ...
            'FontSize', 20, ...
           'MapLatLimit', lat_lim, ...
           'MapLonLimit', lon_lim, ...
           'Frame', 'on', 'Grid', 'on', ...
           'MeridianLabel','on','ParallelLabel','on');
axis off;  

land = shaperead('landareas.shp','UseGeoCoords',true);
geoshow([land.Lat], [land.Lon], ...
        'DisplayType','polygon', ...
        'FaceColor',[0.8 0.8 0.8]);

% 4) Heatmap
surfm(LatG, LonG, ValG,'FaceAlpha',0.5);

for pt = 1:nPts
    textm(lat_vector(pt), lon_vector(pt), ...
          sprintf('%.f', maxGapMinutes(pt)), ...
          'HorizontalAlignment','center', ...
          'VerticalAlignment','middle', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'Color', 'k'); 
end

% 5) Colorbar and Title
cb = colorbar;
caxis([0 60]);
ylabel(cb, 'Max Time (minutes) w/o Coverage');
%title('Heatmap - Lambert Conformal Projection');
set(gca,'FontSize', 25);