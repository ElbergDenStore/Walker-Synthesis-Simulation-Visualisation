function out = p618_aalborg_freq_elevation_study(mode, data_file)
% P618_AALBORG_FREQ_ELEVATION_STUDY Generate and/or plot P.618 attenuation.
%
% This study isolates atmospheric attenuation (At) for one fixed location
% (Aalborg, Denmark) across elevation angle and frequency.
% No FSPL and no orbit-height dependence are included.
%
% Frequencies: 2 GHz, 12 GHz, 20 GHz
% Elevation grid: 5:1:90 deg
% Location: Aalborg (lat=57.0488, lon=9.9217)
%
% Usage:
%   p618_aalborg_freq_elevation_study();                 % generate + plot
%   p618_aalborg_freq_elevation_study("generate");       % generate + save only
%   p618_aalborg_freq_elevation_study("plot");           % plot from saved file only
%   p618_aalborg_freq_elevation_study("plot", "mydata.mat");
%
% Output:
%   Returns a struct with generated/loaded data and figure handle (if plot).

    if nargin < 1 || strlength(string(mode)) == 0
        mode = "both";
    end
    if nargin < 2 || strlength(string(data_file)) == 0
        data_file = "p618_aalborg_freq_elevation.mat";
    end

    mode = lower(string(mode));
    data_file = string(data_file);

    if ~ismember(mode, ["generate", "plot", "both"])
        error("Mode must be one of: 'generate', 'plot', 'both'.");
    end

    settings = struct();
    settings.location_name = "Aalborg, Denmark";
    settings.latitude_deg = 57.0488;
    settings.longitude_deg = 9.9217;
    settings.elevation_deg = 20:1:90;
    settings.frequencies_hz = [2e9, 12e9, 20e9];
    settings.total_annual_exceedance_pct = 1;
    settings.station_height_m = 0;

    D = struct();
    fig = [];

    if mode == "generate" || mode == "both"
        fprintf("\nGenerating P.618 data for %s...\n", settings.location_name);
        D = generate_data(settings);
        save(data_file, "D", "-v7.3");
        fprintf("Saved data: %s\n", data_file);
    end

    if mode == "plot" || mode == "both"
        if isempty(fieldnames(D))
            if ~isfile(data_file)
                error("Data file not found: %s. Run mode='generate' first.", data_file);
            end
            S = load(data_file, "D");
            D = S.D;
            fprintf("Loaded data: %s\n", data_file);
        end
        fig = plot_data(D);
    end

    out = struct("data", D, "figure", fig, "data_file", data_file, "mode", mode);
end

function D = generate_data(settings)
    n_el = numel(settings.elevation_deg);
    n_f = numel(settings.frequencies_hz);

    At_dB = NaN(n_el, n_f);
    Tsky_K = NaN(n_el, n_f);

    total = n_el * n_f;
    done = 0;

    for fi = 1:n_f
        f_hz = settings.frequencies_hz(fi);
        for ei = 1:n_el
            el = settings.elevation_deg(ei);
            [pl, tsky] = safe_p618_point_eval( ...
                settings.latitude_deg, settings.longitude_deg, el, f_hz, ...
                settings.total_annual_exceedance_pct, settings.station_height_m);

            At_dB(ei, fi) = pl.At;
            Tsky_K(ei, fi) = tsky;

            done = done + 1;
            if mod(done, 10) == 0 || done == total
                fprintf("Progress: %6.2f%% (%d/%d)\n", 100 * done / total, done, total);
            end
        end
    end

    D = struct();
    D.location_name = settings.location_name;
    D.latitude_deg = settings.latitude_deg;
    D.longitude_deg = settings.longitude_deg;
    D.elevation_deg = settings.elevation_deg;
    D.frequencies_hz = settings.frequencies_hz;
    D.total_annual_exceedance_pct = settings.total_annual_exceedance_pct;
    D.station_height_m = settings.station_height_m;
    D.At_dB = At_dB;
    D.Tsky_K = Tsky_K;
    D.generated_utc = datetime("now", "TimeZone", "UTC");
end

function fig = plot_data(D)
    fig = figure("Color", "w", "Name", "P.618 Attenuation vs Elevation (Aalborg)");
    hold on;

    colors = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10];

    for fi = 1:numel(D.frequencies_hz)
        plot(D.elevation_deg, D.At_dB(:, fi), "LineWidth", 2.4, "Color", colors(fi, :));
    end

    grid on;
    box on;
    xlabel("Elevation Angle (deg)");
    ylabel("Atmospheric Attenuation At (dB)");

    legend_labels = arrayfun(@(f) sprintf("f = %.0f GHz", f / 1e9), D.frequencies_hz, "UniformOutput", false);
    legend(legend_labels, "Location", "northeast");

    title(sprintf("P.618 Atmospheric Attenuation | %s (%.4f, %.4f) |p=%.2f%%", ...
        D.location_name, D.latitude_deg, D.longitude_deg,D.total_annual_exceedance_pct));

    % subtitle(sprintf("No FSPL, no orbit-height term | Exceedance %.2f%%", D.total_annual_exceedance_pct));
    hold off;

    save_plot_to_screenshots(fig);
end

function save_plot_to_screenshots(fig)
    out_dir = "screenshots";
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    timestamp = datestr(now, "yyyymmdd_HHMMSS");
    base_name = "p618_aalborg_freq_elevation_" + string(timestamp);

    png_path = fullfile(out_dir, base_name + ".png");
    fig_path = fullfile(out_dir, base_name + ".fig");

    exportgraphics(fig, png_path, "Resolution", 300);
    savefig(fig, fig_path);

    fprintf("Saved figure: %s\n", png_path);
end

function [pl, tsky] = safe_p618_point_eval(lat0, lon0, elev, f_hz, exceedance_pct, station_height_m)
    persistent low_freq_note_printed
    if isempty(low_freq_note_printed)
        low_freq_note_printed = false;
    end

    if f_hz < 4e9 && ~low_freq_note_printed
        fprintf("Note: f < 4 GHz triggers an internal P.618 XPD warning; suppressing repeated warning spam during generation.\n");
        low_freq_note_printed = true;
    end

    old_warn = warning('off', 'all');
    cleanup_obj = onCleanup(@() warning(old_warn));

    try
        cfg = p618Config( ...
            "Frequency", f_hz, ...
            "ElevationAngle", elev, ...
            "Latitude", lat0, ...
            "Longitude", lon0, ...
            "TotalAnnualExceedance", exceedance_pct);
        [pl, ~, tsky] = p618PropagationLosses(cfg, "StationHeight", station_height_m);
        return;
    catch ME
        if ~contains(ME.message, "WaterVaporDensity")
            rethrow(ME);
        end
    end

    % Fallback: nearest neighboring points if digital map returns NaN.
    delta_list = [0.25, 0.5, 1, 2, 5];
    lon_dirs = [0, -1, 1, -1, -1, 1, 1, 0, 0];
    lat_dirs = [0, 0, 0, -1, 1, -1, 1, -1, 1];

    for d = delta_list
        for k = 1:numel(lon_dirs)
            lat_try = clamp_lat(lat0 + d * lat_dirs(k));
            lon_try = wrap_lon(lon0 + d * lon_dirs(k));

            cfg = p618Config( ...
                "Frequency", f_hz, ...
                "ElevationAngle", elev, ...
                "Latitude", lat_try, ...
                "Longitude", lon_try, ...
                "TotalAnnualExceedance", exceedance_pct);
            try
                [pl, ~, tsky] = p618PropagationLosses(cfg, "StationHeight", station_height_m);
                return;
            catch ME_try
                if ~contains(ME_try.message, "WaterVaporDensity")
                    rethrow(ME_try);
                end
            end
        end
    end

    error("p618_aalborg_freq_elevation_study:NoValidNeighbor", ...
        "No valid P.618 map value near lat=%.4f lon=%.4f.", lat0, lon0);
end

function lat_c = clamp_lat(lat_in)
    lat_c = max(min(lat_in, 89.999), -89.999);
end

function lon_w = wrap_lon(lon_in)
    lon_w = mod(lon_in + 180, 360) - 180;
end
