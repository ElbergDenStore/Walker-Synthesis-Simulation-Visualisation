function plot_oneweb_tle_altitude()
% PLOT_ONEWEB_TLE_ALTITUDE
%   Propagates 4 OneWeb satellites in the same orbital plane from their
%   TLEs using MATLAB's SGP4 propagator and plots geodetic altitude vs
%   latitude for one orbital pass per satellite.
%
%   All 4 satellites share RAAN ~41.8 deg, inc ~87.89 deg (same plane).
%   TLE epochs: 2026-May-20 UTC
%
%   Usage:
%       plot_oneweb_tle_altitude()

    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'functions'));  % bootstrap so path_setup is found
    path_setup();

    %% TLE definitions ------------------------------------------------------
    sats_def = { ...
        'OW-61609', ...
        '1 61609U 24188R   26140.19895427 -.00000120  00000-0 -37973-3 0  9990', ...
        '2 61609  87.8870  41.8227 0000915 106.0390 254.0838 13.10369549 79625'; ...
        'OW-51629', ...
        '1 51629U 22012H   26140.14084356 -.00000497  00000-0 -14613-2 0  9993', ...
        '2 51629  87.8856  41.8567 0001890  77.6988 282.4351 13.10373994205610'; ...
        'OW-45149', ...
        '1 45149U 20008U   26140.18071808  .00000634  00000-0  17791-2 0  9999', ...
        '2 45149  87.8854  41.8403 0002488  80.3957 279.7451 13.10380156304515'; ...
        'OW-51651', ...
        '1 51651U 22012AF  26140.16244155  .00000302  00000-0  82907-3 0  9990', ...
        '2 51651  87.8855  41.8592 0002020  71.7990 288.3357 13.10375373205654'; ...
    };
    n_sats = size(sats_def, 1);

    %% Use mean motion of first satellite to determine period ---------------
    n_revday = 13.10369549;
    T_s      = 86400 / n_revday;
    fprintf('\nOrbital period: %.2f min\n', T_s/60);

    %% Use latest TLE epoch as scenario start (minimises age-of-data) ------
    epoch_strs = {'26140.19895427','26140.14084356','26140.18071808','26140.16244155'};
    epoch_dts  = cellfun(@(s) parse_tle_epoch(s), epoch_strs);
    start_dt   = max(epoch_dts);
    fprintf('Scenario start: %s UTC\n\n', char(start_dt, 'yyyy-MM-dd HH:mm:ss'));

    %% Write all 4 TLEs to one temp file ------------------------------------
    tmp_tle = [tempname, '.tle'];
    fid = fopen(tmp_tle, 'w');
    for k = 1:n_sats
        fprintf(fid, '%s\n%s\n%s\n', sats_def{k,1}, sats_def{k,2}, sats_def{k,3});
    end
    fclose(fid);

    %% Build scenario -------------------------------------------------------
    sc = satelliteScenario;
    sc.StartTime  = start_dt;
    sc.StopTime   = start_dt + seconds(T_s + 30);
    sc.SampleTime = 5;

    sat_objs = satellite(sc, tmp_tle, 'OrbitPropagator','sgp4');
    delete(tmp_tle);

    %% Colours (distinct, colourblind-friendly) -----------------------------
    cols = [0.00 0.45 0.74;   % blue
            0.85 0.33 0.10;   % orange
            0.47 0.67 0.19;   % green
            0.49 0.18 0.56];  % purple

    %% Plot -----------------------------------------------------------------
    f = figure('Color','w', 'Position',[100 100 820 460]);
    hold on;

    for k = 1:n_sats
        P   = squeeze(states(sat_objs(k), 'CoordinateFrame','ECEF'))';  % Nx3 m
        lla = ecef2lla(P);
        lat = lla(:,1);
        alt = lla(:,3) / 1e3;

        % Ascending pass, sorted by latitude
        asc_mask    = [false; diff(lat) >= 0];
        [lat_s, ix] = sort(lat(asc_mask));
        alt_tmp     = alt(asc_mask);
        alt_s       = alt_tmp(ix);

        fprintf('%s  min=%.2f km  max=%.2f km  p-p=%.1f m\n', ...
            sats_def{k,1}, min(alt), max(alt), (max(alt)-min(alt))*1e3);

        plot(lat_s, alt_s, '-', 'Color', cols(k,:), 'LineWidth', 1.8, ...
            'DisplayName', sats_def{k,1});
    end

    xlabel('Geodetic latitude (deg)');
    ylabel('Geodetic altitude (km)');
    title('OneWeb – same-plane satellites, altitude vs latitude (SGP4, 2026-May-20)', ...
        'FontWeight','bold');
    legend('Location','best');
    grid on;
    xlim([-90 90]);

    %% Save -----------------------------------------------------------------
    script_dir = fileparts(mfilename('fullpath'));
    out_dir    = fullfile(script_dir, '..', 'figures');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    fname = fullfile(out_dir, 'oneweb_same_plane_altitude_vs_latitude.png');
    exportgraphics(f, fname, 'Resolution',300);
    fprintf('\nPlot saved -> %s\n', fname);
end

%--------------------------------------------------------------------------
function dt = parse_tle_epoch(epoch_str)
% Parse TLE epoch string 'YYDDD.FFFFFFFF' into a datetime (UTC)
    yy        = str2double(epoch_str(1:2));
    year_full = 2000 + yy;
    day_frac  = str2double(epoch_str(3:end));
    dt        = datetime(year_full, 1, 1, 'TimeZone','UTC') + days(day_frac - 1);
end
