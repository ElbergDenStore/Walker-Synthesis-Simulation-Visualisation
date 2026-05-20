function constellation_altitude_spread(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG)
% CONSTELLATION_ALTITUDE_SPREAD
%   Propagates every satellite in an asymmetrical Walker Star with BOTH the
%   SGP4 and the numerical propagator and shows the min/max geodetic altitude
%   across ALL satellites as a function of latitude for each.
%   The spread reveals how bad the initialisation error (e=0, w=0 passed as
%   osculating elements) is and whether it differs between propagators.
%
%   Usage:
%       constellation_altitude_spread()                          % defaults
%       constellation_altitude_spread(1000, 87, 65, 5, 2.5)

    if nargin < 1 || isempty(HEIGHT_KM),      HEIGHT_KM      = 1000; end
    if nargin < 2 || isempty(INC_DEG),        INC_DEG        = 90;   end
    if nargin < 3 || isempty(PLANES),         PLANES         = 2;   end
    if nargin < 4 || isempty(SATS_PER_PLANE), SATS_PER_PLANE = 5;    end
    if nargin < 5 || isempty(MIN_EL_DEG),     MIN_EL_DEG     = 20;  end

    total_sats = PLANES * SATS_PER_PLANE;
    fprintf('\n=== Constellation altitude spread analysis ===\n');
    fprintf('  h = %g km   i = %g deg\n', HEIGHT_KM, INC_DEG);
    fprintf('  Planes = %d   Sats/plane = %d   Total = %d\n', ...
        PLANES, SATS_PER_PLANE, total_sats);
    fprintf('  Seam min-elevation = %g deg\n\n', MIN_EL_DEG);

    %% Orbit period ----------------------------------------------------------
    Re_m = 6378.14e3;
    mu   = 3.986004418e14;
    a_m  = Re_m + HEIGHT_KM*1e3;
    T_s  = 2*pi * sqrt(a_m^3 / mu);
    fprintf('  Orbital period: %.2f min\n\n', T_s/60);

    %% Run propagators -------------------------------------------------------
    [bins_kep,  lat_ctrs, diag_kep ] = run_propagator(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG, T_s, 'two-body-keplerian', 'Kepler');
    [bins_sgp4, ~,        diag_sgp4] = run_propagator(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG, T_s, 'sgp4',               'SGP4');
    [bins_num,  ~,        diag_num ] = run_propagator(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG, T_s, 'numerical',          'numerical');
    [bins_rk4,  ~,        diag_rk4 ] = run_rk4_propagator(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG, T_s);

    %% Side-by-side comparison table ----------------------------------------
    print_comparison_table({'Kepler','SGP4','numerical','RK4+J2 (own)'}, {diag_kep, diag_sgp4, diag_num, diag_rk4});

    script_dir = fileparts(mfilename('fullpath'));
    out_dir    = fullfile(script_dir, 'figures');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    %% Plot Kepler -----------------------------------------------------------
    plot_and_save(bins_kep, lat_ctrs, total_sats, HEIGHT_KM, INC_DEG, ...
        'Kepler (two-body)', [0.49 0.18 0.56], out_dir, ...
        sprintf('altitude_spread_kepler_%dp%ds_%dkm_i%.0f.png', PLANES, SATS_PER_PLANE, HEIGHT_KM, INC_DEG));

    %% Plot SGP4 -------------------------------------------------------------
    plot_and_save(bins_sgp4, lat_ctrs, total_sats, HEIGHT_KM, INC_DEG, ...
        'SGP4', [0.85 0.33 0.10], out_dir, ...
        sprintf('altitude_spread_sgp4_%dp%ds_%dkm_i%.0f.png', PLANES, SATS_PER_PLANE, HEIGHT_KM, INC_DEG));

    %% Plot numerical --------------------------------------------------------
    plot_and_save(bins_num, lat_ctrs, total_sats, HEIGHT_KM, INC_DEG, ...
        'numerical', [0.00 0.45 0.74], out_dir, ...
        sprintf('altitude_spread_numerical_%dp%ds_%dkm_i%.0f.png', PLANES, SATS_PER_PLANE, HEIGHT_KM, INC_DEG));

    %% Plot RK4 + J2 (own propagator) ---------------------------------------
    plot_and_save(bins_rk4, lat_ctrs, total_sats, HEIGHT_KM, INC_DEG, ...
        'RK4+J2 (own)', [0.47 0.67 0.19], out_dir, ...
        sprintf('altitude_spread_rk4j2_%dp%ds_%dkm_i%.0f.png', PLANES, SATS_PER_PLANE, HEIGHT_KM, INC_DEG));
end

%% =========================================================================
function plot_and_save(bins, lat_ctrs, total_sats, HEIGHT_KM, INC_DEG, label, col, out_dir, fname)
    f = figure('Color','w','Position',[80 80 980 520]);

    % Top panel: min/max/median altitude vs latitude
    ax1 = subplot(2,1,1);
    hold on;
    fill([lat_ctrs, fliplr(lat_ctrs)], ...
         [bins.alt_max, fliplr(bins.alt_min)], ...
         col + (1-col)*0.55, 'EdgeColor','none', 'DisplayName','min–max range');
    plot(lat_ctrs, bins.alt_med, 'k-', 'LineWidth',1.2, 'DisplayName','median');
    yline(HEIGHT_KM, ':r', 'LineWidth',1.0, 'DisplayName','requested h');
    xlabel('Geodetic latitude (deg)');
    ylabel('Geodetic altitude (km)');
    title(sprintf('%s altitude spread across %d-sat constellation  (h=%g km, i=%g deg)', ...
        label, total_sats, HEIGHT_KM, INC_DEG), 'FontWeight','bold');
    legend('Location','best');
    grid on; xlim([-90 90]);

    % Bottom panel: peak-to-peak spread per latitude bin
    ax2 = subplot(2,1,2);
    bar(lat_ctrs, bins.spread_m, 1, 'FaceColor',col, 'EdgeColor','none');
    xlabel('Geodetic latitude (deg)');
    ylabel('Altitude spread (m)');
    title('Peak-to-peak altitude spread within each 1-degree latitude bin');
    grid on; xlim([-90 90]);

    linkaxes([ax1, ax2], 'x');

    exportgraphics(f, fullfile(out_dir, fname), 'Resolution',200);
    fprintf('  Plot saved -> %s\n', fullfile(out_dir, fname));
end

%% =========================================================================
function [bins, lat_ctrs, diag] = run_propagator(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG, T_s, prop, label)
% Builds the constellation with the given propagator, propagates one orbit,
% bins altitude by 1-degree latitude, and returns summary statistics.

    fprintf('  [%s] Building + propagating %d satellites...\n', ...
        label, PLANES*SATS_PER_PLANE);

    sc = satelliteScenario;
    sc.StartTime  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
    sc.StopTime   = sc.StartTime + seconds(T_s + 30);
    sc.SampleTime = 30;

    tic;
    sats = asymmetrical_walker_star_generation(sc, HEIGHT_KM*1e3, INC_DEG, ...
        PLANES, SATS_PER_PLANE, MIN_EL_DEG, prop);
    t_build = toc;
    fprintf('  [%s] Satellite generation:  %.2f s\n', label, t_build);

    tic;
    Nsat   = numel(sats);
    sat_lat_cell = cell(Nsat,1);
    sat_alt_cell = cell(Nsat,1);
    for k = 1:Nsat
        P   = squeeze(states(sats(k), 'CoordinateFrame','ECEF'))';
        lla = ecef2lla(P);
        sat_lat_cell{k} = lla(:,1);
        sat_alt_cell{k} = lla(:,3)/1e3;
    end
    all_lat = vertcat(sat_lat_cell{:});
    all_alt = vertcat(sat_alt_cell{:});
    t_prop  = toc;
    fprintf('  [%s] Propagation + ECEF->LLA: %.2f s\n', label, t_prop);

    [bins, lat_ctrs] = bin_altitudes(all_lat, all_alt);
    diag = compute_diagnostics(sat_lat_cell, sat_alt_cell, HEIGHT_KM, label);
end

%% =========================================================================
function [bins, lat_ctrs, diag] = run_rk4_propagator(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG, T_s)
% Same output schema as run_propagator but uses the in-house RK4+J2
% integrator (1 s step) in ECI, with eci2ecef + ecef2lla for the geometry.

    label = 'RK4+J2';
    fprintf('  [%s] Building + propagating %d satellites (dt = 1 s)...\n', ...
        label, PLANES*SATS_PER_PLANE);

    start_time  = datetime('1-Jun-2025 12:00:00','TimeZone','UTC');
    duration_s  = T_s + 30;
    dt          = 1.0;
    sample_step = 30.0;

    [all_lat, all_alt, t_build, t_prop] = propagate_walker_star_rk4_j2( ...
        HEIGHT_KM*1e3, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG, ...
        start_time, duration_s, dt, sample_step);

    fprintf('  [%s] Satellite generation:  %.2f s\n', label, t_build);
    fprintf('  [%s] Propagation + ECI->ECEF->LLA: %.2f s\n', label, t_prop);

    % Reshape flat arrays back into per-satellite cells.
    % Layout in the RK4 helper is: rows = (sIdx-1)*Nsat + (1:Nsat), so each
    % block of Nsat rows is ONE time sample across all sats.
    Nsat       = PLANES * SATS_PER_PLANE;
    n_samples  = numel(all_lat) / Nsat;
    lat_mat    = reshape(all_lat, Nsat, n_samples);   % Nsat x n_samples
    alt_mat    = reshape(all_alt, Nsat, n_samples);
    sat_lat_cell = cell(Nsat,1);
    sat_alt_cell = cell(Nsat,1);
    for k = 1:Nsat
        sat_lat_cell{k} = lat_mat(k,:).';
        sat_alt_cell{k} = alt_mat(k,:).';
    end

    [bins, lat_ctrs] = bin_altitudes(all_lat, all_alt);
    diag = compute_diagnostics(sat_lat_cell, sat_alt_cell, HEIGHT_KM, label);
end

%% =========================================================================
function [bins, lat_ctrs] = bin_altitudes(all_lat, all_alt)
    lat_edges = -90:1:90;
    lat_ctrs  = lat_edges(1:end-1) + 0.5;
    n_bins    = numel(lat_ctrs);

    alt_min = NaN(1, n_bins);
    alt_max = NaN(1, n_bins);
    alt_med = NaN(1, n_bins);
    for b = 1:n_bins
        mask = all_lat >= lat_edges(b) & all_lat < lat_edges(b+1);
        if any(mask)
            v = all_alt(mask);
            alt_min(b) = min(v);
            alt_max(b) = max(v);
            alt_med(b) = median(v);
        end
    end

    bins.alt_min  = alt_min;
    bins.alt_max  = alt_max;
    bins.alt_med  = alt_med;
    bins.spread_m = (alt_max - alt_min) * 1e3;
end

%% =========================================================================
function diag = compute_diagnostics(sat_lat_cell, sat_alt_cell, HEIGHT_KM, label)
% Per-satellite diagnostics:
%   - altitude at the ascending equator crossing
%   - altitude at the descending equator crossing
%   - altitude at the north pole pass (lat closest to +90)
%   - altitude at the south pole pass (lat closest to -90)
% From these we derive:
%   - asc-desc asymmetry  (per sat)
%   - N-S asymmetry       (per sat)
%   - inter-sat spread    (across constellation at each "probe" latitude)

    Nsat = numel(sat_lat_cell);
    alt_asc   = NaN(Nsat,1);
    alt_desc  = NaN(Nsat,1);
    alt_north = NaN(Nsat,1);
    alt_south = NaN(Nsat,1);
    perigee   = NaN(Nsat,1);
    apogee    = NaN(Nsat,1);

    for k = 1:Nsat
        lat = sat_lat_cell{k};
        alt = sat_alt_cell{k};
        if numel(lat) < 3, continue; end

        % Sign-change in latitude => equator crossings.
        s    = sign(lat);
        cross_idx = find(s(1:end-1) .* s(2:end) < 0);

        asc_alts  = [];
        desc_alts = [];
        for ci = cross_idx.'
            % Linear interp on time between samples ci and ci+1 to lat = 0.
            l1 = lat(ci); l2 = lat(ci+1);
            a1 = alt(ci); a2 = alt(ci+1);
            f  = l1 / (l1 - l2);          % fraction to zero crossing
            a0 = a1 + f * (a2 - a1);
            if l2 > l1                     % ascending (going N)
                asc_alts(end+1) = a0;  %#ok<AGROW>
            else
                desc_alts(end+1) = a0; %#ok<AGROW>
            end
        end
        if ~isempty(asc_alts),  alt_asc(k)  = mean(asc_alts);  end
        if ~isempty(desc_alts), alt_desc(k) = mean(desc_alts); end

        % North / south pole approach (highest |lat| reached).
        [~, in_]  = max(lat);
        [~, is_]  = min(lat);
        alt_north(k) = alt(in_);
        alt_south(k) = alt(is_);

        perigee(k) = min(alt);
        apogee(k)  = max(alt);
    end

    diag.label       = label;
    diag.alt_asc     = alt_asc;
    diag.alt_desc    = alt_desc;
    diag.alt_north   = alt_north;
    diag.alt_south   = alt_south;
    diag.perigee     = perigee;
    diag.apogee      = apogee;
    diag.HEIGHT_KM   = HEIGHT_KM;

    asc_desc_m = (alt_asc  - alt_desc ) * 1e3;
    ns_m       = (alt_north - alt_south) * 1e3;
    pa_m       = (apogee   - perigee  ) * 1e3;

    diag.asc_desc_m = asc_desc_m;
    diag.ns_m       = ns_m;
    diag.pa_m       = pa_m;

    fprintf('\n  ---- Diagnostics [%s] ----\n', label);
    fprintf('  Requested altitude              : %8.3f km\n', HEIGHT_KM);
    fprintf('  Per-sat asc-desc asymmetry  (m) : mean %+8.2f   std %7.2f   |max| %7.2f\n', ...
        mean(asc_desc_m,'omitnan'), std(asc_desc_m,'omitnan'), max(abs(asc_desc_m),[],'omitnan'));
    fprintf('  Per-sat N-S asymmetry       (m) : mean %+8.2f   std %7.2f   |max| %7.2f\n', ...
        mean(ns_m,'omitnan'),       std(ns_m,'omitnan'),       max(abs(ns_m),[],'omitnan'));
    fprintf('  Per-sat apogee-perigee      (m) : mean %8.2f   std %7.2f   max  %7.2f\n', ...
        mean(pa_m,'omitnan'),       std(pa_m,'omitnan'),       max(pa_m,[],'omitnan'));
    fprintf('  -- inter-sat spread at same latitude (km vs requested h) --\n');
    fprintf('  Asc  equator: min %+7.3f   max %+7.3f   spread %7.2f m\n', ...
        min(alt_asc)-HEIGHT_KM,   max(alt_asc)-HEIGHT_KM,   (max(alt_asc)-min(alt_asc))*1e3);
    fprintf('  Desc equator: min %+7.3f   max %+7.3f   spread %7.2f m\n', ...
        min(alt_desc)-HEIGHT_KM,  max(alt_desc)-HEIGHT_KM,  (max(alt_desc)-min(alt_desc))*1e3);
    fprintf('  North pass  : min %+7.3f   max %+7.3f   spread %7.2f m\n', ...
        min(alt_north)-HEIGHT_KM, max(alt_north)-HEIGHT_KM, (max(alt_north)-min(alt_north))*1e3);
    fprintf('  South pass  : min %+7.3f   max %+7.3f   spread %7.2f m\n', ...
        min(alt_south)-HEIGHT_KM, max(alt_south)-HEIGHT_KM, (max(alt_south)-min(alt_south))*1e3);
end

%% =========================================================================
function print_comparison_table(labels, diags)
    fprintf('\n========== PROPAGATOR COMPARISON SUMMARY ==========\n');
    fprintf('%-15s | %12s | %12s | %12s | %14s\n', ...
        'Propagator', 'asc-desc(m)', 'N-S(m)', 'apo-peri(m)', 'eq-spread(m)');
    fprintf('%s\n', repmat('-',1,75));
    for k = 1:numel(diags)
        d = diags{k};
        fprintf('%-15s | %12.2f | %12.2f | %12.2f | %14.2f\n', ...
            labels{k}, ...
            max(abs(d.asc_desc_m),[],'omitnan'), ...
            max(abs(d.ns_m),[],'omitnan'),       ...
            max(d.pa_m,[],'omitnan'),            ...
            (max(d.alt_asc)-min(d.alt_asc))*1e3);
    end
    fprintf('%s\n', repmat('=',1,75));
    fprintf('Expectation: all columns ~0 if the constellation is initialised correctly.\n');
    fprintf('Columns are max absolute values across all satellites.\n\n');
end

