function max_P_tx_dBm = PFD_calc(target_PFD_MHz, g_tx, bandwidth, orbit_height, min_elevation_UE)
% PFD_CALC  Maximum transmit power allowed by a Power Flux Density limit.
%   Back-solves the per-beam transmit power MAX_P_TX_DBM (dBm) that just meets the
%   PFD limit TARGET_PFD_MHz at the worst-case nadir slant range, given transmit
%   gain G_TX (dBi), BANDWIDTH (Hz), ORBIT_HEIGHT (m) and the minimum user
%   elevation MIN_ELEVATION_UE (deg).
    Re = 6371e3;    
    cos_exponent = 1.5;

    % Calculate Slant Range 
    slant_range = -Re.*sind(min_elevation_UE) + sqrt(Re^2.*sind(min_elevation_UE).^2 - (Re^2-(Re+orbit_height)^2));
    
    % Calculate Steering Angle and Loss
    sin_theta_tx = (Re * cosd(min_elevation_UE)) / (Re+orbit_height);
    theta_tx = asind(sin_theta_tx);
    loss_steer = -10 * log10(cosd(theta_tx)^cos_exponent);
    
    % Calculate Area
    area_m2 = 4 * pi * (slant_range).^2;
    
    % Calculate tx power
    % PFD_1MHz = g_tx + P_tx_dBm - loss_steer - 10*log10(area_m2) - 10*log10(B/1e6);
    max_P_tx_dBm = target_PFD_MHz + 30 - g_tx + loss_steer + 10*log10(area_m2) + 10*log10(bandwidth/1e6);
end