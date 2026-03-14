function G_tx = get_adjusted_tx_gain(orbit_height, min_elevation_UE, frequency)
    %%% To achieve same size beams on earth for different orbit heights, tx
    %%% gain can be adjusted

    % A reference is needed
    f_ref = 20e9;   
    h_ref = 1200e3; 
    G_tx_ref = 42;  
    Re = 6371e3;    
    
    % Calculate Slant Range for the Reference Altitude
    slant_ref = -Re*sind(min_elevation_UE) + sqrt(Re^2*sind(min_elevation_UE)^2 - (Re^2-(Re+h_ref)^2));
    
    % Calculate Slant Range for the CURRENT Altitude in the sweep
    current_h = orbit_height;
    slant_current = -Re*sind(min_elevation_UE) + sqrt(Re^2*sind(min_elevation_UE)^2 - (Re^2-(Re+current_h)^2));
    
    % Scale the Antenna Gain to keep the Ground Footprint constant
    G_tx = G_tx_ref + 20 * log10(slant_current / slant_ref) - 10*log10(f_ref/frequency);
end