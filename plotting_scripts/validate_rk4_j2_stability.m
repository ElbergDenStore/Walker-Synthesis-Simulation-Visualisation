function validate_rk4_j2_stability(HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG)
% VALIDATE_RK4_J2_STABILITY
%   Propagates a Walker Star constellation for 250 simulation hours using
%   the standalone RK4+J2 propagator and checks for orbital energy drift and
%   inter-satellite spread growth over time.
%
%   Metrics checked
%     1. Geocentric radius trend  - running mean of |r(t)| for satellite 1.
%        Should be flat; any slope indicates energy drift (integrator failure).
%     2. Ascending equator altitude stability - altitude at every equator
%        crossing (ascending pass) for all satellites over 250 h.
%        All points should lie within a narrow band (< 50 m) for all time.
%     3. Inter-satellite equator spread per orbit - max minus min altitude
%        at ascending crossings within each orbital period.  Should stay
%        near the initial value (~10 m) and not grow.
%     4. Geocentric radius oscillation - |r(t)| for satellite 1 over the
%        first 3 orbital periods, confirming the frozen orbit shape.
%
%   Usage
%     validate_rk4_j2_stability()
%     validate_rk4_j2_stability(1000, 90, 5, 13, 20)

    if nargin < 1 || isempty(HEIGHT_KM),      HEIGHT_KM      = 1000; end
    if nargin < 2 || isempty(INC_DEG),        INC_DEG        = 90;   end
    if nargin < 3 || isempty(PLANES),         PLANES         = 5;    end
    if nargin < 4 || isempty(SATS_PER_PLANE), SATS_PER_PLANE = 13;   end
    if nargin < 5 || isempty(MIN_EL_DEG),     MIN_EL_DEG     = 20;   end

    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

    DURATION_H  = 250;
    dt          = 1.0;
    sample_step = 660.0;
    duration_s  = DURATION_H * 3600;
    start_time  = datetime('1-Jun-2025 12:00:00', 'TimeZone', 'UTC');

    Re_m  = 6378.137e3;
    mu    = 3.986004418e14;
    a_m   = Re_m + HEIGHT_KM * 1e3;
    T_orb = 2*pi * sqrt(a_m^3 / mu);   % orbital period (s)
    Nsat  = PLANES * SATS_PER_PLANE;

    fprintf('\n=== RK4+J2 250-hour stability validation ===\n');
    fprintf('  h = %g km   i = %g deg   %d × %d constellation (%d sats)\n', ...
        HEIGHT_KM, INC_DEG, PLANES, SATS_PER_PLANE, Nsat);
    fprintf('  Duration : %d h  |  Orbital periods : %.1f\n\n', ...
        DURATION_H, duration_s / T_orb);

    %% --- 1. Generate frozen-orbit initial states ---
    fprintf('  Generating frozen-orbit initial states ...\n');
    tic;
    [r0, v0] = generate_walker_star_states(HEIGHT_KM*1e3, INC_DEG, PLANES, SATS_PER_PLANE, MIN_EL_DEG);
    fprintf('  IC generation : %.2f s\n\n', toc);

    %% --- 2. Propagate ---
    fprintf('  Propagating %d satellites × %d hours (dt = %g s, sample = %g s) ...\n', ...
        Nsat, DURATION_H, dt, sample_step);
    tic;
    [r_ecef_hist, ~] = propagate_rk4_j2(r0, v0, start_time, duration_s, dt, sample_step);
    t_prop = toc;
    n_samples = size(r_ecef_hist, 3);
    fprintf('  Propagation   : %.1f s  (%d samples)\n\n', t_prop, n_samples);

    %% --- 3. Convert all ECEF to LLA in one vectorised call ---
    fprintf('  Converting ECEF → LLA ...\n');
    tic;
    % Layout: permute to [Nsat × n_samples × 3], reshape to [(Nsat*n_samples) × 3]
    r_flat  = reshape(permute(r_ecef_hist, [1, 3, 2]), Nsat*n_samples, 3);
    lla_flat = ecef2lla(r_flat);
    % lat_mat(k,i) = geodetic latitude of satellite k at sample i
    lat_mat  = reshape(lla_flat(:,1), Nsat, n_samples);
    alt_mat  = reshape(lla_flat(:,3), Nsat, n_samples) / 1e3;   % km
    fprintf('  LLA conversion : %.1f s\n\n', toc);

    t_h = (0:n_samples-1) * sample_step / 3600;   % time axis in hours

    %% --- 4. Geocentric radius for satellite 1 ---
    % |r_ecef| = |r_eci| because ECEF is a pure rotation of ECI.
    r_mag_km = squeeze(sqrt(sum(r_ecef_hist(1,:,:).^2, 2))) / 1e3;   % n_samples × 1

    % Running mean over one orbital period — reveals any secular drift
    win      = max(1, round(T_orb / sample_step));
    r_mean   = movmean(r_mag_km, win);

    %% --- 5. Find ascending equator crossings for all satellites ---
    fprintf('  Extracting ascending equator crossings ...\n');
    asc_t_h  = [];   % time of crossing (h)
    asc_alt  = [];   % geodetic altitude (km) at crossing
    asc_sat  = [];   % satellite index

    for k = 1:Nsat
        lat_k = lat_mat(k,:);
        alt_k = alt_mat(k,:);
        s = sign(lat_k);
        for i = 1:n_samples-1
            if s(i) < 0 && s(i+1) > 0   % ascending: south → north
                f = lat_k(i) / (lat_k(i) - lat_k(i+1));
                asc_t_h(end+1)  = t_h(i)  + f * (t_h(i+1)  - t_h(i));  %#ok<AGROW>
                asc_alt(end+1)  = alt_k(i) + f * (alt_k(i+1) - alt_k(i));  %#ok<AGROW>
                asc_sat(end+1)  = k;  %#ok<AGROW>
            end
        end
    end

    %% --- 6. Per-orbit inter-satellite spread at ascending equator ---
    orbit_edges   = 0 : T_orb/3600 : DURATION_H;
    n_orbits      = numel(orbit_edges) - 1;
    spread_m      = NaN(n_orbits, 1);
    orbit_t_mid   = (orbit_edges(1:end-1) + orbit_edges(2:end)) / 2;

    for oi = 1:n_orbits
        mask = asc_t_h >= orbit_edges(oi) & asc_t_h < orbit_edges(oi+1);
        if sum(mask) > 1
            spread_m(oi) = (max(asc_alt(mask)) - min(asc_alt(mask))) * 1e3;
        end
    end

    %% --- 7. Print summary ---
    early  = asc_t_h < 10;
    late   = asc_t_h > (DURATION_H - 10);
    spread_early = (max(asc_alt(early)) - min(asc_alt(early))) * 1e3;
    spread_late  = (max(asc_alt(late))  - min(asc_alt(late)))  * 1e3;
    r_drift_m    = (r_mean(end) - r_mean(1)) * 1e3;    % metres over 250 h

    fprintf('\n  ── Stability summary ─────────────────────────────────────────\n');
    fprintf('  Metric                              Value\n');
    fprintf('  ─────────────────────────────────────────────────────────────\n');
    fprintf('  Geocentric radius drift (250 h)     %+.2f m\n', r_drift_m);
    fprintf('  Asc-equator spread, t = 0–10 h      %.1f m\n', spread_early);
    fprintf('  Asc-equator spread, t = 240–250 h   %.1f m\n', spread_late);
    fprintf('  Spread change (late – early)         %+.1f m\n', spread_late - spread_early);
    fprintf('  Median per-orbit spread              %.1f m\n', median(spread_m,'omitnan'));
    fprintf('  Max    per-orbit spread              %.1f m\n', max(spread_m,[],'omitnan'));
    fprintf('  ─────────────────────────────────────────────────────────────\n\n');

    %% --- 8. Plots ---
    script_dir = fileparts(mfilename('fullpath'));
    out_dir    = fullfile(script_dir, 'figures');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    fig = figure('Color','w', 'Position', [60 60 1200 820]);

    %  Panel 1 (top-left): geocentric radius over 250 h  ──────────────────
    ax1 = subplot(2, 2, 1);
    plot(t_h, r_mag_km, 'Color', [0.6 0.6 0.8], 'LineWidth', 0.5); hold on;
    plot(t_h, r_mean,   'b-', 'LineWidth', 1.5);
    xlabel('Time (h)');
    ylabel('Geocentric radius (km)');
    title('Satellite 1 — geocentric radius');
    legend({'instantaneous', sprintf('running mean (T_{orb} window)')}, 'Location', 'best');
    grid on;

    %  Panel 2 (top-right): altitude at ascending equator, all sats  ───────
    ax2 = subplot(2, 2, 2);
    cmap = parula(Nsat);
    for k = 1:Nsat
        mask_k = asc_sat == k;
        scatter(asc_t_h(mask_k), (asc_alt(mask_k) - HEIGHT_KM)*1e3, 4, ...
            cmap(k,:), 'filled', 'MarkerFaceAlpha', 0.6);
        hold on;
    end
    xlabel('Time (h)');
    ylabel('Altitude offset from nominal (m)');
    title(sprintf('Ascending equator altitude — all %d sats', Nsat));
    yline(0, ':k', 'LineWidth', 1);
    grid on;
    colormap(ax2, parula(Nsat));
    cb = colorbar(ax2);
    cb.Label.String = 'Satellite index';
    clim([1, Nsat]);

    %  Panel 3 (bottom-left): per-orbit spread  ───────────────────────────
    ax3 = subplot(2, 2, 3);
    bar(orbit_t_mid, spread_m, 1, 'FaceColor', [0.47 0.67 0.19], 'EdgeColor', 'none');
    xlabel('Time (h)');
    ylabel('Inter-sat spread (m)');
    title('Ascending equator spread per orbital period');
    grid on;
    ylim([0, max(max(spread_m,[],'omitnan')*1.5, 50)]);

    %  Panel 4 (bottom-right): first 3 orbits — orbital shape  ────────────
    ax4 = subplot(2, 2, 4);
    n_3orb = min(n_samples, round(3 * T_orb / sample_step) + 1);
    plot(t_h(1:n_3orb), r_mag_km(1:n_3orb) - a_m/1e3, ...
        'b-', 'LineWidth', 1.2);
    xlabel('Time (h)');
    ylabel('r – a_{nominal} (km)');
    title('Satellite 1 — radial excursion, first 3 orbits');
    grid on;

    sgtitle(sprintf('RK4+J2 stability validation — %d×%d Walker Star, h=%g km, i=%g°, %d h', ...
        PLANES, SATS_PER_PLANE, HEIGHT_KM, INC_DEG, DURATION_H), ...
        'FontWeight', 'bold', 'FontSize', 12);

    fname = sprintf('stability_rk4j2_%dp%ds_%dkm_i%.0f_%dh.png', ...
        PLANES, SATS_PER_PLANE, HEIGHT_KM, INC_DEG, DURATION_H);
    exportgraphics(fig, fullfile(out_dir, fname), 'Resolution', 200);
    fprintf('  Figure saved → %s\n', fullfile(out_dir, fname));
end
