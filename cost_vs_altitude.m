clear all
close all
clc

%% Sweep parameter

min_elevation_angle_degree = [20];
min_elevation_angle = deg2rad(min_elevation_angle_degree);
sat_altitude_km = [550:50:1200];

r_km = zeros(length(sat_altitude_km), length(min_elevation_angle));

for i = 1:length(sat_altitude_km)
    for j = 1:length(min_elevation_angle)
        
        r_km(i,j) = cono_cobertura_function(sat_altitude_km(i),min_elevation_angle(j));

    end
end

%% N Satelites

A_km = pi*r_km.^2;

N_satellites = ceil(510072000./A_km); %Earth Surface = 510 072 000 km²​​​ 
N_satellites = N_satellites./54; %Normalization assuming we need 54 satellites at 1000 km

%% Estimated Cost
                          % 550  600  650  700  750  800  850  900  950  1000 1050 1100 1150 1200
normalized_cost_550_1200 = [1.10 1.20 1.30 1.40 1.50 1.62 1.75 1.90 2.05 2.21 2.38 2.57 2.77 3.00]./2.21; %Normalized to 1000 km


normalized_satellites_multiplied_cost = normalized_cost_550_1200.*N_satellites';

figure
set(gca,'FontSize',25)
grid on
grid minor
yyaxis left
plot(sat_altitude_km, normalized_cost_550_1200, 'LineWidth', 3, 'LineStyle','-', 'Marker','o');
xlim([550, 1200])
ylim([0.4, 1.4])
xlabel('Satellite Altitude (km)')
ylabel('Normalized Cost /Satellite')
yyaxis right
plot(sat_altitude_km, normalized_satellites_multiplied_cost(1,:), 'LineWidth', 3, 'LineStyle','-', 'Marker','o');
xlim([550, 1200])
ylim([0.625, 1.25])
xlabel('Satellite Altitude (km)')
ylabel('Normalized Total Cost')
