function LUT = generate_p618_lookup_table(totalAnnualExceedance, frequency_hz, lat_limits, lon_limits)
% GENERATE_P618_LOOKUP_TABLE Precompute regional/global P.618 attenuation LUT.
%
% Builds a 3D lookup with:
%   latitude  : lat_limits(1):1:lat_limits(2) (deg)
%   longitude : lon_limits(1):1:lon_limits(2) (deg)
%   elevation : 10:10:90   (deg)
%
% Outputs:
%   LUT.At_dB(lat, lon, el)  : Total atmospheric attenuation [dB]
%   LUT.Tsky_K(lat, lon, el) : Sky noise temperature [K]
%
% The output is saved to a MAT file for very fast interpolation use.
%
% Example:
%   generate_p618_lookup_table();
%   generate_p618_lookup_table(1);
%   generate_p618_lookup_table(1, 12e9);
%   generate_p618_lookup_table(1, 13.5e9, [54.58, 83.67], [-73.17, 33.50]);

    if nargin < 1
        totalAnnualExceedance = 1;
    end
    if nargin < 2
        frequency_hz = 13.5e9;
    end
    if nargin < 3
        lat_limits = [-90, 90];
    end
    if nargin < 4
        lon_limits = [-180, 180];
    end

    validateattributes(lat_limits, {'double'}, {'vector', 'numel', 2, 'real', 'finite'});
    validateattributes(lon_limits, {'double'}, {'vector', 'numel', 2, 'real', 'finite'});
    lat_limits = sort(lat_limits(:)');
    lon_limits = sort(lon_limits(:)');

    if lat_limits(1) < -90 || lat_limits(2) > 90
        error("generate_p618_lookup_table:InvalidLatLimits", "lat_limits must be within [-90, 90].");
    end
    if lon_limits(1) < -180 || lon_limits(2) > 180
        error("generate_p618_lookup_table:InvalidLonLimits", "lon_limits must be within [-180, 180].");
    end

    lat_grid = lat_limits(1):1:lat_limits(2);
    lon_grid = lon_limits(1):1:lon_limits(2);
    el_grid = 10:10:90;

    n_lat = numel(lat_grid);
    n_lon = numel(lon_grid);
    n_el = numel(el_grid);
    n_points = n_lat * n_lon;

    is_global = isequal(lat_limits, [-90, 90]) && isequal(lon_limits, [-180, 180]);
    if is_global
        fprintf("\nGenerating global P.618 lookup table...\n");
    else
        fprintf("\nGenerating regional P.618 lookup table...\n");
    end
    fprintf("Lat range: [%.3f, %.3f] deg | Lon range: [%.3f, %.3f] deg\n", ...
        lat_limits(1), lat_limits(2), lon_limits(1), lon_limits(2));
    fprintf("Grid: %d lat x %d lon x %d el = %d cells\n", n_lat, n_lon, n_el, n_lat * n_lon * n_el);

    At_dB = NaN(n_lat, n_lon, n_el, "single");
    Tsky_K = NaN(n_lat, n_lon, n_el, "single");

    % if isempty(gcp("nocreate"))
    %     parpool("fullspeed");
    % end

    t_start = tic;
    dq = parallel.pool.DataQueue;
    points_done = 0;
    afterEach(dq, @on_progress_points);
    fprintf("Progress:   0.00%% (0/%d lat-lon points) | ETA: --\n", n_points);

    parfor i_lat = 1:n_lat
        warning("off", "all");

        % PRE-ALLOCATE: Create ONE config object per latitude slice (worker)
        cfg = p618Config("Frequency", frequency_hz, "TotalAnnualExceedance", totalAnnualExceedance);

        lat = lat_grid(i_lat);
        At_slice = zeros(n_lon, n_el, "single");
        Tsky_slice = zeros(n_lon, n_el, "single");
        lon_update_step = max(1, floor(n_lon / 20));
        last_reported_lon = 0;

        for i_lon = 1:n_lon
            lon = lon_grid(i_lon);
            for i_el = 1:n_el
                el = el_grid(i_el);
                
                % Pass the pre-allocated object into the eval function
                [pl, tsky] = safe_p618_eval_fast(cfg, lat, lon, el);
                
                At_slice(i_lon, i_el) = single(pl.At);
                Tsky_slice(i_lon, i_el) = single(tsky);
            end

            if mod(i_lon, lon_update_step) == 0 || i_lon == n_lon
                delta_points = i_lon - last_reported_lon;
                send(dq, delta_points);
                last_reported_lon = i_lon;
            end
        end

        At_dB(i_lat, :, :) = At_slice;
        Tsky_K(i_lat, :, :) = Tsky_slice;
    end

    LUT = struct();
    LUT.lat_deg = lat_grid;
    LUT.lon_deg = lon_grid;
    LUT.el_deg = el_grid;
    LUT.frequency_hz = frequency_hz;
    LUT.total_annual_exceedance_pct = totalAnnualExceedance;
    LUT.station_height_m = 0;
    LUT.lat_limits_deg = lat_limits;
    LUT.lon_limits_deg = lon_limits;
    LUT.At_dB = At_dB;
    LUT.Tsky_K = Tsky_K;
    LUT.generated_utc = datetime("now", "TimeZone", "UTC");
    LUT.grid_shape = [n_lat, n_lon, n_el];

    out_name = sprintf("p618_lookup_%0.1f.mat", frequency_hz / 1e9);
    save(out_name, "LUT", "-v7.3");

    fprintf("Done in %.1f s\n", toc(t_start));
    fprintf("Saved: %s\n", out_name);

    function on_progress_points(delta_points)
        points_done = points_done + delta_points;
        elapsed = toc(t_start);
        frac = points_done / n_points;

        if frac > 0
            eta_sec = elapsed * (1 / frac - 1);
            eta_str = format_eta(eta_sec);
        else
            eta_str = "--";
        end

        fprintf("Progress: %6.2f%% (%d/%d lat-lon points) | ETA: %s\n", ...
            100 * frac, points_done, n_points, eta_str);
        drawnow limitrate;
    end
end

function eta_str = format_eta(seconds_left)
    if ~isfinite(seconds_left) || seconds_left < 0
        eta_str = "--";
        return;
    end

    h = floor(seconds_left / 3600);
    m = floor(mod(seconds_left, 3600) / 60);
    s = floor(mod(seconds_left, 60));
    eta_str = sprintf('%02dh:%02dm:%02ds', h, m, s);
end

function [pl, tsky] = safe_p618_eval_fast(cfg, lat0, lon0, elev)
    % Update the properties of the existing object
    cfg.Latitude = lat0;
    cfg.Longitude = lon0;
    cfg.ElevationAngle = elev;

    try
        [pl, ~, tsky] = p618PropagationLosses(cfg, "StationHeight", 0);
        return;
    catch ME
        if ~contains(ME.message, "WaterVaporDensity")
            rethrow(ME);
        end
    end

    % Fallback: retry nearest neighboring coordinates when digital maps are
    % undefined at the exact point.
    delta_list = [0.25, 0.5, 1, 2, 5];
    lon_dirs = [0, -1, 1, -1, -1, 1, 1, 0, 0];
    lat_dirs = [0, 0, 0, -1, 1, -1, 1, -1, 1];

    for d = delta_list
        for k = 1:numel(lon_dirs)
            lat_try = clamp_lat(lat0 + d * lat_dirs(k));
            lon_try = wrap_lon(lon0 + d * lon_dirs(k));

            % Update coordinates for the fallback attempt
            cfg.Latitude = lat_try;
            cfg.Longitude = lon_try;

            try
                [pl, ~, tsky] = p618PropagationLosses(cfg, "StationHeight", 0);
                return;
            catch ME_try
                if ~contains(ME_try.message, "WaterVaporDensity")
                    rethrow(ME_try);
                end
            end
        end
    end

    error("generate_p618_lookup_table:NoValidNeighbor", ...
        "No valid P.618 map value near lat=%.3f lon=%.3f el=%.1f deg.", lat0, lon0, elev);
end

function lat_c = clamp_lat(lat_in)
    lat_c = max(min(lat_in, 89.999), -89.999);
end

function lon_w = wrap_lon(lon_in)
    lon_w = mod(lon_in + 180, 360) - 180;
end
