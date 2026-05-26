Re = 6378.14e3;     
cos_exponent = 1.5;
el_vec = linspace(20,90,100);
orbit_height = [1000e3];

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'plotting_scripts/figures', 'steering_loss');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

for i = 1:length(orbit_height)
    R_sat = Re + orbit_height(i);
    theta_rx = 90 - el_vec; 
    sin_theta_tx = (Re .* cosd(el_vec)) ./ R_sat;
    theta_tx = asind(sin_theta_tx); 
    
    
    loss_tx = -10 * log10(cosd(theta_tx).^cos_exponent);

    
    loss_rx = -10 * log10(cosd(theta_rx).^cos_exponent);


    f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 500 400]);

    set(0, 'DefaultAxesFontSize', 14);
    plot(el_vec,loss_tx, "DisplayName", 'Satellite')
    hold on;
    grid on;
    plot(el_vec,loss_rx, "DisplayName", 'UE')
    plot(el_vec,loss_rx+loss_tx,"--", "DisplayName", 'total')
    legend('Location', 'northeast');
    xlabel("Elevation (deg)")
    ylabel("Steering Loss (dB)")
    title(sprintf("Steering loss @ %d km, for cos()^{1.5}",orbit_height(i)/1e3))

    out_file = fullfile(out_dir, sprintf('steering_loss_%dkm.png', round(orbit_height(i)/1e3)));
    exportgraphics(f1, out_file, 'Resolution', 300);
end