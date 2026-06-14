% orbit_height_km = 550;
% target_radius_km = 16;
% gain_from_nadir_beam_size(orbit_height_km, target_radius_km)

function G_tx = gain_from_nadir_beam_size(orbit_height_km, target_radius_km)
    % gain_from_nadir_beam_size Determines the required antenna boresight 
    % gain to achieve a specific footprint radius on Earth.
    %
    % Inputs:
    %   orbit_height_km  - The altitude of the satellite in kilometers
    %   target_radius_km - The desired radius of the beam footprint on the ground
    %
    % Outputs:
    %   G_tx             - The required antenna gain in dBi

    % --- 1. Constants ---
    Re = 6371; % Earth's mean radius in km

    % --- 2. Pure Geometric Derivation ---
    % Calculate the Earth Central Angle (gamma) in radians
    gamma_rad = target_radius_km / Re;
    
    % Calculate the distance from Earth's center to the satellite
    r_sat = Re + orbit_height_km;
    
    % Calculate Slant Range (d) to the edge of the footprint via Law of Cosines
    slant_range = sqrt(Re^2 + r_sat^2 - 2 * Re * r_sat * cos(gamma_rad));
    
    % Calculate the Nadir Half-Beamwidth (eta) via Law of Sines
    eta_rad = asin((Re * sin(gamma_rad)) / slant_range);
    
    % The full beamwidth is simply twice the half-beamwidth
    beamwidth_rad = 2 * eta_rad;
    
    % --- 3. Calculate Required Gain ---
    % Use the standard directivity approximation for an ideal aperture (G = 4*pi / theta^2)
    G_tx_linear = (4 * pi) / (beamwidth_rad^2);
    G_tx = 10 * log10(G_tx_linear); % Convert linear gain to dBi
end