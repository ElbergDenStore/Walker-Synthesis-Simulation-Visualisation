% amount of beams from  
% function BeamGrid = calculate_hexagonal_beams(G_tx_dBi, f_Hz, orbit_height_m, Min_Elev_deg, target_EIRP_or_PFD, FRF)

% =========================================================================
% SATELLITE & RF GEOMETRY UTILITY FORMULAS
% =========================================================================
repo_root = fileparts(mfilename('fullpath'));
while ~isfile(fullfile(repo_root, 'functions', 'path_setup.m')), repo_root = fileparts(repo_root); end
addpath(fullfile(repo_root, 'functions'));
path_setup();
freq_hz = 12e9;
efficiency = 0.6;
element_gain_dbi = 6;

orbit_height_km = 1200;
target_radius_km = 30;
gain_dbi = gain_from_nadir_beam_size(orbit_height_km, target_radius_km)


r = calc_coverage_radius(550, 20);
hpbw_deg = calc_hpbw_from_gain(gain_dbi)
size = calc_antenna_size_from_gain(gain_dbi, freq_hz, efficiency)
elements = calc_array_elements(gain_dbi, element_gain_dbi)


orbit_height_km = 8000;
% target_radius_km = 16;
gain_dbi = gain_from_nadir_beam_size(orbit_height_km, target_radius_km)

r = calc_coverage_radius(8000, 20);
hpbw_deg = calc_hpbw_from_gain(gain_dbi)
size = calc_antenna_size_from_gain(gain_dbi, freq_hz, efficiency)
elements = calc_array_elements(gain_dbi, element_gain_dbi)
sub_arrays = elements / 16;



% -------------------------------------------------------------------------
% 1. Coverage Radius from Altitude and Min Elevation
% -------------------------------------------------------------------------
function r_km = calc_coverage_radius(h_km, min_elev_deg)
    Re = 6371; % Earth radius in km
    elev_rad = deg2rad(min_elev_deg);
    
    % Law of Sines to find the Nadir Angle (eta)
    eta_rad = asin((Re / (Re + h_km)) * cos(elev_rad));
    
    % Earth Central Angle (gamma)
    gamma_rad = (pi/2) - elev_rad - eta_rad;
    
    % Arc length on the surface of the Earth
    r_km = Re * gamma_rad;
end

% -------------------------------------------------------------------------
% 2. Coverage Area from Altitude and Min Elevation
% -------------------------------------------------------------------------
function area_km2 = calc_coverage_area(h_km, min_elev_deg)
    Re = 6371; 
    elev_rad = deg2rad(min_elev_deg);
    
    eta_rad = asin((Re / (Re + h_km)) * cos(elev_rad));
    gamma_rad = (pi/2) - elev_rad - eta_rad;
    
    % Area of a spherical cap
    area_km2 = 2 * pi * Re^2 * (1 - cos(gamma_rad));
end

% -------------------------------------------------------------------------
% 3. HPBW Approximation from Antenna Gain
% -------------------------------------------------------------------------
function hpbw_deg = calc_hpbw_from_gain(gain_dbi)
    % Assuming an ideal uniform circular aperture
    g_linear = 10^(gain_dbi / 10);
    hpbw_rad = sqrt((4 * pi) / g_linear);
    
    hpbw_deg = rad2deg(hpbw_rad);
end

% -------------------------------------------------------------------------
% 4. Antenna Size (Diameter) from Gain and Frequency
% -------------------------------------------------------------------------
function diameter_m = calc_antenna_size_from_gain(gain_dbi, freq_hz, efficiency)
    if nargin < 3
        efficiency = 0.6; % Default to typical parabolic reflector efficiency
    end
    
    c = 299792458; % Speed of light in m/s
    lambda = c / freq_hz;
    g_linear = 10^(gain_dbi / 10);
    
    % Derived from G = efficiency * (pi * D / lambda)^2
    diameter_m = (lambda / pi) * sqrt(g_linear / efficiency);
end

% -------------------------------------------------------------------------
% 5. Beam Footprint from Gain, Altitude, and Elevation Angle
% -------------------------------------------------------------------------
function [area_km2, semi_major_km, semi_minor_km] = calc_beam_footprint(gain_dbi, h_km, elev_deg)
    % Calculates the elliptical footprint on the ground when the beam 
    % is steered to a specific elevation angle.
    
    Re = 6371;
    elev_rad = deg2rad(elev_deg);
    
    % Get HPBW in radians
    g_linear = 10^(gain_dbi / 10);
    hpbw_rad = sqrt((4 * pi) / g_linear);
    
    % Calculate slant range to the target elevation via Law of Cosines
    slant_km = -Re * sin(elev_rad) + sqrt((Re * sin(elev_rad))^2 + 2*Re*h_km + h_km^2);
    
    % Calculate minor axis (across-track width)
    semi_minor_km = (slant_km * hpbw_rad) / 2;
    
    % Calculate major axis (along-track width, stretched by the grazing angle)
    semi_major_km = semi_minor_km / sin(elev_rad);
    
    % Area of the resulting ellipse
    area_km2 = pi * semi_major_km * semi_minor_km;
end

% -------------------------------------------------------------------------
% 6. Approximate Array Elements from Antenna Gain
% -------------------------------------------------------------------------
function num_elements = calc_array_elements(gain_dbi, element_gain_dbi)
    if nargin < 2
        element_gain_dbi = 5; % Typical directivity for a patch antenna
    end
    
    % Array Gain = Element Gain + 10*log10(N)
    % Solved for N:
    num_elements = ceil(10^((gain_dbi - element_gain_dbi) / 10));
end