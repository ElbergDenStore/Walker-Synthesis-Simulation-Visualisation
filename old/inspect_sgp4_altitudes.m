% inspect_sgp4_altitudes.m
% Compares SGP4 altitude variation for Walker Star vs Walker Delta.
%
% The geodetic offset vs nominal is geometric: orbits use a spherical Earth
% (r=6378.14 km) but ecef2lla uses WGS84. At high inclinations satellites
% spend time near the poles where WGS84 surface is ~21 km closer to Earth's
% centre, so geodetic altitude reads higher than nominal.
% Plot 4 shows orbital radius |r| which removes the WGS84 geodetic effect.

clear; close all; clc;

%% ── 1. Run both constellations ──────────────────────────────────────────
HEIGHT_KM = 1200;
DURATION  = "long";   % "short" ~1 orbit; "medium" ~12 h; "long" ~48 h

fprintf('=== Walker Star ===\n');
Cfg_star  = get_cfg(HEIGHT_KM, "walkerstar",  "small", DURATION);
R_star    = analyse_constellation(Cfg_star);

fprintf('\n=== Walker Delta ===\n');
Cfg_delta = get_cfg(HEIGHT_KM, "walkerdelta", "small", DURATION);
R_delta   = analyse_constellation(Cfg_delta);

%% ── 2. Comparison plots ─────────────────────────────────────────────────

%% Plot 1: Fleet altitude envelopes side by side
figure('Name','Altitude Envelopes','Position',[50 100 1300 420]);
for k = 1:2
    if k==1; R=R_star;  label='Walker Star';  col=[0.4 0.6 0.9];
    else;    R=R_delta; label='Walker Delta'; col=[0.9 0.6 0.4]; end
    subplot(1,2,k);
    fill([R.time_min(:); flipud(R.time_min(:))], ...
         [R.min_alt_t(:); flipud(R.max_alt_t(:))], col, ...
         'EdgeColor','none','DisplayName','Fleet min/max');
    hold on;
    plot(R.time_min, R.mean_alt_t, 'Color', col*0.7, 'LineWidth', 1.5, 'DisplayName','Fleet mean');
    yline(HEIGHT_KM,'r--','LineWidth',1.2,'DisplayName',sprintf('Nominal %d km',HEIGHT_KM));
    xlabel('Time (min)'); ylabel('Geodetic altitude (km)');
    title(sprintf('%s  [%d sats, %d km nominal]', label, R.num_sats, HEIGHT_KM));
    legend('Location','best'); grid on;
end

%% Plot 2: Per-satellite mean altitude deviation
figure('Name','Mean Altitude Deviation','Position',[50 560 1300 360]);
for k = 1:2
    if k==1; R=R_star;  label='Walker Star';  col=[0.4 0.6 0.9];
    else;    R=R_delta; label='Walker Delta'; col=[0.9 0.6 0.4]; end
    subplot(1,2,k);
    bar(sort(R.mean_alt_per_sat - HEIGHT_KM), 'FaceColor', col);
    hold on; yline(0,'r--','LineWidth',1.2); grid on;
    xlabel('Satellite (sorted)'); ylabel('Mean alt deviation (km)');
    title(sprintf('%s — offset: %.3f km, spread: %.4f km', ...
        label, mean(R.mean_alt_per_sat)-HEIGHT_KM, range(R.mean_alt_per_sat)));
end

%% Plot 3: Per-satellite peak-to-peak oscillation
figure('Name','Peak-to-Peak Oscillation','Position',[50 100 1300 360]);
for k = 1:2
    if k==1; R=R_star;  label='Walker Star';  col=[0.4 0.6 0.9];
    else;    R=R_delta; label='Walker Delta'; col=[0.9 0.6 0.4]; end
    subplot(1,2,k);
    bar(sort(R.pp_alt_per_sat), 'FaceColor', col);
    xlabel('Satellite (sorted)'); ylabel('Peak-to-peak (km)');
    title(sprintf('%s — mean p-p: %.3f km', label, mean(R.pp_alt_per_sat)));
    grid on;
end

%% Plot 4: Orbital radius |r| - r_earth  (removes WGS84 geodetic bias)
r_earth_km = 6378.14;
figure('Name','Orbital Radius (no WGS84 bias)','Position',[50 560 1300 420]);
for k = 1:2
    if k==1; R=R_star;  label='Walker Star';  col=[0.4 0.6 0.9];
    else;    R=R_delta; label='Walker Delta'; col=[0.9 0.6 0.4]; end
    radius_alt = R.orbital_radius_km - r_earth_km;  % [nT x nSats]
    min_r = min(radius_alt,[],2); max_r = max(radius_alt,[],2); mean_r = mean(radius_alt,2);
    subplot(1,2,k);
    fill([R.time_min(:); flipud(R.time_min(:))], ...
         [min_r(:); flipud(max_r(:))], col, 'EdgeColor','none','DisplayName','min/max');
    hold on;
    plot(R.time_min, mean_r, 'Color', col*0.7, 'LineWidth', 1.5, 'DisplayName','mean');
    yline(HEIGHT_KM,'r--','LineWidth',1.2,'DisplayName','Nominal |r|');
    xlabel('Time (min)'); ylabel('|r| - r_{earth} (km)');
    title(sprintf('%s — orbital radius alt (p-p: %.4f km)', label, range(radius_alt(:))));
    legend('Location','best'); grid on;
end

%% ── Helper function ─────────────────────────────────────────────────────
function R = analyse_constellation(Cfg)
    r_earth = 6378.14e3;
    fprintf('Constellation: %d planes × %d sats = %d total | inc=%.1f°\n', ...
        Cfg.Num_planes, Cfg.Sats_per_plane, Cfg.Total_sats, Cfg.Inclination);

    sc = satelliteScenario;
    sc.StartTime  = Cfg.StartTime;
    sc.StopTime   = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;

    if Cfg.WalkerStar
        sats = asymmetrical_walker_star_generation(sc, Cfg.Orbit_height, ...
                   Cfg.Inclination, Cfg.Num_planes, Cfg.Sats_per_plane, ...
                   Cfg.Min_elevation_UE);
    else
        sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
                   Cfg.Inclination, Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
                   Name="InspectSats", OrbitPropagator="sgp4");
    end

    fprintf('Propagating %d satellites...\n', numel(sats));
    [sat_pos_raw, ~, simTimes] = states(sats, "CoordinateFrame", "ECEF");
    % sat_pos_raw: [3 x nT x nSats]

    num_sats = numel(sats);
    nT       = size(sat_pos_raw, 2);

    % Geodetic altitude via ecef2lla [nT x nSats]
    xyz_flat = reshape(permute(sat_pos_raw, [2 3 1]), num_sats * nT, 3);
    lla_flat = ecef2lla(xyz_flat);
    alt_km   = reshape(lla_flat(:,3)/1e3, nT, num_sats);

    % Orbital radius |r| — no WGS84 geodetic effect [nT x nSats]
    xyz_norm          = sqrt(sum(reshape(sat_pos_raw, 3, []).^2, 1));
    orbital_radius_km = reshape(xyz_norm/1e3, nT, num_sats);

    mean_alt_per_sat = mean(alt_km, 1);
    pp_alt_per_sat   = max(alt_km,[],1) - min(alt_km,[],1);
    nominal          = Cfg.Orbit_height/1e3;

    fprintf('  Nominal (spherical)      : %.3f km\n', nominal);
    fprintf('  Geodetic mean (fleet)    : %.3f km  (σ = %.4f km)\n', mean(mean_alt_per_sat), std(mean_alt_per_sat));
    fprintf('  Spread in mean alt       : range = %.4f km\n', range(mean_alt_per_sat));
    fprintf('  Mean peak-to-peak (geo)  : %.3f km\n', mean(pp_alt_per_sat));
    r_alt = orbital_radius_km - r_earth/1e3;
    fprintf('  Orbital radius alt mean  : %.3f km  (p-p across fleet: %.4f km)\n', mean(r_alt(:)), range(r_alt(:)));

    R.num_sats          = num_sats;
    R.time_min          = minutes(simTimes - simTimes(1));
    R.alt_km            = alt_km;
    R.orbital_radius_km = orbital_radius_km;
    R.mean_alt_per_sat  = mean_alt_per_sat;
    R.pp_alt_per_sat    = pp_alt_per_sat;
    R.min_alt_t         = min(alt_km,[],2);
    R.max_alt_t         = max(alt_km,[],2);
    R.mean_alt_t        = mean(alt_km,2);
end
