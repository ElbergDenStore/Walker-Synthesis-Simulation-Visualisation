function out = benchmark_aer_methods(height_km, ue_grid_size, duration)
% BENCHMARK_AER_METHODS Compare MATLAB AER workflow against vectorized geometry math.
%
% This benchmark runs both methods on the same constellation and UEs,
% then prints timing and numerical agreement. No plotting is performed.
%
% Example:
%   benchmark_aer_methods();
%   benchmark_aer_methods(1000, "medium", "short");

    if nargin < 1
        height_km = 1000;
    end
    if nargin < 2
        ue_grid_size = "big";
    end
    if nargin < 3
        duration = "long";
    end

    % Ensure repo root + all functions/ subfolders are on the path,
    % regardless of the caller's current working directory.
    repo_root = fileparts(fileparts(mfilename('fullpath')));  % validators/ -> repo root
    addpath(fullfile(repo_root, 'functions'));                % bootstrap so path_setup is found
    path_setup();

    Cfg = get_cfg(height_km, "walkerdelta", ue_grid_size, duration);

    fprintf('\n=== AER vs Vectorized Geometry Benchmark ===\n');
    fprintf('Height: %d km | UEs: %d | Sats: %d\n', ...
        height_km, numel(Cfg.Flat_UE_array.Lats), Cfg.Total_sats);

    UE_lats = Cfg.Flat_UE_array.Lats(:);
    UE_lons = Cfg.Flat_UE_array.Lons(:);
    num_ues = numel(UE_lats);

    %% Method 1: MATLAB aer()
    fprintf('\nRunning aer() method...\n');
    tic;
    [num_visible_aer, best_range_aer, best_el_aer, best_az_aer, best_sat_aer, aer_setup_time, aer_prop_time] = ...
        run_aer_benchmark(Cfg, UE_lats, UE_lons);
    t_aer = toc;

    %% Method 2: Vectorized ENU math
    fprintf('Running vectorized method...\n');
    tic;
    [num_visible_vec, best_range_vec, best_el_vec, best_az_vec, best_sat_vec, vec_setup_time, vec_prop_time] = ...
        run_vectorized_benchmark(Cfg, UE_lats, UE_lons);
    t_vec = toc;

    %% Compare results
    [max_range_err, max_el_err, max_az_err, sat_id_mismatch_pct, vis_mismatch_pct] = ...
        compare_results(best_range_aer, best_range_vec, best_el_aer, best_el_vec, ...
                        best_az_aer, best_az_vec, best_sat_aer, best_sat_vec, ...
                        num_visible_aer, num_visible_vec);

    speedup = t_aer / t_vec;

    fprintf('\n--- Timing ---\n');
    fprintf('aer() end-to-end time    : %.3f s\n', t_aer);
    fprintf('Vectorized end-to-end    : %.3f s\n', t_vec);
    fprintf('Speedup (aer/vectorized) : %.2fx\n', speedup);
    fprintf('AER setup time          : %.3f s\n', aer_setup_time);
    fprintf('AER geometry time       : %.3f s\n', aer_prop_time);
    fprintf('Vector setup/state time : %.3f s\n', vec_setup_time);
    fprintf('Vector geometry time    : %.3f s\n', vec_prop_time);

    fprintf('\n--- Numerical Agreement ---\n');
    fprintf('Max |Range error| (m)     : %.6g\n', max_range_err);
    fprintf('Max |Elevation error| (deg): %.6g\n', max_el_err);
    fprintf('Max |Azimuth error| (deg) : %.6g\n', max_az_err);
    fprintf('Best Sat ID mismatch (%%)  : %.6f\n', sat_id_mismatch_pct);
    fprintf('Visible-count mismatch (%%): %.6f\n', vis_mismatch_pct);

    out = struct;
    out.config = struct("height_km", height_km, "ue_grid_size", ue_grid_size, "duration", duration);
    out.size = struct("numUEs", num_ues, "numSats", Cfg.Total_sats);
    out.timing = struct( ...
        "aer_seconds", t_aer, ...
        "vectorized_seconds", t_vec, ...
        "speedup", speedup, ...
        "aer_setup_seconds", aer_setup_time, ...
        "aer_propagation_seconds", aer_prop_time, ...
        "vector_setup_seconds", vec_setup_time, ...
        "vector_propagation_seconds", vec_prop_time);
    out.accuracy = struct( ...
        "max_range_error_m", max_range_err, ...
        "max_elevation_error_deg", max_el_err, ...
        "max_azimuth_error_deg", max_az_err, ...
        "sat_id_mismatch_percent", sat_id_mismatch_pct, ...
        "visible_count_mismatch_percent", vis_mismatch_pct);
end

function [num_visible, best_range, best_el, best_az, best_sat, setup_time, propagation_time] = run_aer_benchmark(Cfg, UE_lats, UE_lons)
    setup_tic = tic;
    sc = satelliteScenario;
    sc.StartTime = Cfg.StartTime;
    sc.StopTime = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;

    r_earth = 6378.137e3;
    sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");

    ue = groundStation(sc, UE_lats, UE_lons);
    setup_time = toc(setup_tic);

    prop_tic = tic;
    [num_visible, best_range, best_el, best_az, best_sat] = run_aer_method(ue, sats, Cfg.Min_elevation_UE);
    propagation_time = toc(prop_tic);
end

function [num_visible, best_range, best_el, best_az, best_sat] = run_aer_method(ue, sats, min_elev)
    num_ues = numel(ue);

    % Get dimensions from first UE to avoid expensive dynamic growth.
    [~, el0, ~] = aer(ue(1), sats);
    nT = size(el0, 2);

    num_visible = zeros(num_ues, nT);
    best_range = NaN(num_ues, nT);
    best_el = NaN(num_ues, nT);
    best_az = NaN(num_ues, nT);
    best_sat = NaN(num_ues, nT);

    for idx = 1:num_ues
        [az_mat, el_mat, r_mat] = aer(ue(idx), sats);

        valid_mask = el_mat >= min_elev;
        num_visible(idx, :) = sum(valid_mask, 1);
        has_service = any(valid_mask, 1);

        r_temp = r_mat;
        r_temp(~valid_mask) = Inf;
        [best_ranges, best_sat_idx] = min(r_temp, [], 1);

        if any(has_service)
            best_range(idx, has_service) = best_ranges(has_service);
            best_sat(idx, has_service) = best_sat_idx(has_service);

            valid_cols = find(has_service);
            num_rows = size(el_mat, 1);
            best_sats_valid = best_sat_idx(has_service);
            lin_idxs = best_sats_valid + (valid_cols - 1) * num_rows;

            best_el(idx, has_service) = el_mat(lin_idxs);
            best_az(idx, has_service) = az_mat(lin_idxs);
        end
    end
end

function [num_visible, best_range, best_el, best_az, best_sat, setup_time, propagation_time] = run_vectorized_benchmark(Cfg, UE_lats, UE_lons)
    setup_tic = tic;
    sc = satelliteScenario;
    sc.StartTime = Cfg.StartTime;
    sc.StopTime = Cfg.StopTime;
    sc.SampleTime = Cfg.SampleTime;

    r_earth = 6378.137e3;
    sats = walkerDelta(sc, Cfg.Orbit_height + r_earth, ...
        Cfg.Inclination, Cfg.Total_sats, Cfg.Num_planes, Cfg.Phasing, ...
        Name="S4D", OrbitPropagator="sgp4");

    [sat_pos_raw, ~, ~] = states(sats, "CoordinateFrame", "ECEF");
    sat_pos_ecef = permute(sat_pos_raw, [1, 3, 2]);
    setup_time = toc(setup_tic);

    prop_tic = tic;
    [num_visible, best_range, best_el, best_az, best_sat] = run_vectorized_method(UE_lats, UE_lons, sat_pos_ecef, Cfg.Min_elevation_UE);
    propagation_time = toc(prop_tic);
end

function [num_visible, best_range, best_el, best_az, best_sat] = run_vectorized_method(UE_lats, UE_lons, sat_pos_ecef, min_elev)
    ue_pos_ecef = lla2ecef([UE_lats, UE_lons, zeros(numel(UE_lats), 1)]);

    num_ues = numel(UE_lats);
    num_sats = size(sat_pos_ecef, 2);
    nT = size(sat_pos_ecef, 3);

    num_visible = zeros(num_ues, nT);
    best_range = NaN(num_ues, nT);
    best_el = NaN(num_ues, nT);
    best_az = NaN(num_ues, nT);
    best_sat = NaN(num_ues, nT);

    for idx = 1:num_ues
        ue_xyz = ue_pos_ecef(idx, :)';
        lat = UE_lats(idx);
        lon = UE_lons(idx);

        vec_ecef = sat_pos_ecef - ue_xyz;

        slat = sind(lat);
        clat = cosd(lat);
        slon = sind(lon);
        clon = cosd(lon);

        R_ecef_to_enu = [
            -slon,           clon,          0;
            -slat*clon,     -slat*slon,     clat;
             clat*clon,      clat*slon,     slat
        ];

        vec_ecef_flat = reshape(vec_ecef, 3, []);
        vec_enu_flat = R_ecef_to_enu * vec_ecef_flat;
        vec_enu = reshape(vec_enu_flat, 3, num_sats, nT);

        E = reshape(vec_enu(1, :, :), num_sats, nT);
        N = reshape(vec_enu(2, :, :), num_sats, nT);
        U = reshape(vec_enu(3, :, :), num_sats, nT);

        r_mat = sqrt(E.^2 + N.^2 + U.^2);
        el_mat = asind(U ./ r_mat);
        az_mat = atan2d(E, N);
        az_mat(az_mat < 0) = az_mat(az_mat < 0) + 360;

        valid_mask = el_mat >= min_elev;
        num_visible(idx, :) = sum(valid_mask, 1);
        has_service = any(valid_mask, 1);

        r_temp = r_mat;
        r_temp(~valid_mask) = Inf;
        [best_ranges, best_sat_idx] = min(r_temp, [], 1);

        if any(has_service)
            best_range(idx, has_service) = best_ranges(has_service);
            best_sat(idx, has_service) = best_sat_idx(has_service);

            valid_cols = find(has_service);
            best_sats_valid = best_sat_idx(has_service);
            lin_idxs = best_sats_valid + (valid_cols - 1) * num_sats;

            best_el(idx, has_service) = el_mat(lin_idxs);
            best_az(idx, has_service) = az_mat(lin_idxs);
        end
    end
end

function [max_range_err, max_el_err, max_az_err, sat_id_mismatch_pct, vis_mismatch_pct] = ...
    compare_results(best_range_aer, best_range_vec, best_el_aer, best_el_vec, ...
                    best_az_aer, best_az_vec, best_sat_aer, best_sat_vec, ...
                    num_visible_aer, num_visible_vec)

    common_valid = ~isnan(best_range_aer) & ~isnan(best_range_vec);
    if any(common_valid(:))
        max_range_err = max(abs(best_range_aer(common_valid) - best_range_vec(common_valid)));
        max_el_err = max(abs(best_el_aer(common_valid) - best_el_vec(common_valid)));

        d_az = abs(best_az_aer(common_valid) - best_az_vec(common_valid));
        d_az = min(d_az, 360 - d_az);
        max_az_err = max(d_az);
    else
        max_range_err = NaN;
        max_el_err = NaN;
        max_az_err = NaN;
    end

    same_service_mask = ~isnan(best_sat_aer) & ~isnan(best_sat_vec);
    if any(same_service_mask(:))
        sat_id_mismatch_pct = 100 * mean(best_sat_aer(same_service_mask) ~= best_sat_vec(same_service_mask));
    else
        sat_id_mismatch_pct = NaN;
    end

    vis_mismatch_pct = 100 * mean(num_visible_aer(:) ~= num_visible_vec(:));
end
