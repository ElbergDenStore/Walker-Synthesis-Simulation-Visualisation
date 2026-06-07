function BeamGrid = calculate_oneweb_beams()
% CALCULATE_ONEWEB_BEAMS  Build the OneWeb phased-array beam grid (16 fan beams).
%   Returns a BeamGrid struct describing the OneWeb 2x32-element phased array:
%   element/array gains, the cosine steering-loss exponent, and the 16 beam
%   centre directions spaced so their 3 dB points just touch.  Takes no inputs.
    Re_km = 6378.137; % Earth radius in km, WGS84
    orbit_height_km = 1200; 

    BeamGrid.OneWeb = true; % Bit quirky loss calculation to the mechanical presteering.

    % OneWeb Phased Array specs
    BeamGrid.Nu = 2;  % 2 element in u-direction (wide beam, e.g., cross-track). Only 1 row can be seen, but the beam specs in FCC filing corresponds to 2 spaced 0.5 lambda based on the contour plots
    BeamGrid.Nv = 32; % 32 elements in v-direction (narrow fan beam, e.g., along-track) 
    beam_count = 16;
    BeamGrid.num_beams = beam_count;

    BeamGrid.Cos_exponent = 1.5; % loss when steering
    BeamGrid.Element_gain = 6; %dBi
    BeamGrid.Max_gain = 10*log10(BeamGrid.Nu * BeamGrid.Nv) + BeamGrid.Element_gain; % Max gain for plotting and power derivations

    % Calculate the exact 3dB beamwidth in v-space
    % 1.391 is the constant for the 3dB point of a sinc function
    v_3dB_width = 2 * (1.391 / (BeamGrid.Nv * (pi/2))); % ≈ 0.0553
    
    % Determine steering centers (v_center)
    % Space the 16 beams so their 3dB points touch exactly
    max_v = ((beam_count - 1) / 2) * v_3dB_width; % ≈ 0.4147
    
    BeamGrid.v_center = linspace(-max_v, max_v, beam_count);
    BeamGrid.u_center = zeros(1, beam_count); 
    
    % Sanity Check: asind(max_v) will be approximately 24.5 degrees
    
    % Calculate Slant Range to Earth for each beam center
    theta_rad = asin(BeamGrid.v_center);
    R_sat = Re_km + orbit_height_km;
    
    slant_range_m = (R_sat .* cos(theta_rad) - sqrt(Re_km^2 - (R_sat .* sin(theta_rad)).^2)) .* 1000;
    
    % Calculate EIRP per beam to meet FCC PFD limits
    target_pfd_dBmHz = -140.5 + 30 - 10*log10(4000); % ≈ -146.52 dBm/Hz/m^2
    area_spreading_dB = 10 * log10(4 * pi * (slant_range_m.^2));
    
    BeamGrid.BeamCenter_EIRP_dBmHz = target_pfd_dBmHz + area_spreading_dB;
    
    % Build Neighbor Matrix 
    BeamGrid.neighbor_idx = nan(beam_count, 2);
    for i = 1:beam_count
        neighbor_list = [];
        if i > 1
            neighbor_list(end+1) = i - 1;
        end
        if i < beam_count
            neighbor_list(end+1) = i + 1;
        end
        
        pad_size = 2 - length(neighbor_list);
        BeamGrid.neighbor_idx(i, :) = [neighbor_list, nan(1, pad_size)];
    end
end