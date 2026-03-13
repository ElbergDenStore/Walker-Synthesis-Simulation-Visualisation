function [r_km] = cono_cobertura_function(sat_altitude_km,min_elevation_angle_rad)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sat_altitude = sat_altitude_km*1000;

Re = 6378135; % Earth's radius in meters
Rs = Re + sat_altitude;

a = Rs + Re;
b = 2 * Rs / tan(pi/2 - min_elevation_angle_rad);
c = -sat_altitude;

D = b^2 - 4*a*c;
t1 = (-b + sqrt(D)) / (2 * a);

angle = 2 * atan(t1); % Angle in radians
    
% Compute surface coverage radius
r = Re * angle;

r_km = r/1000;
end