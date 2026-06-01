function figure_handles = plot_beamgrid(BeamGrid)
%PLOT_BEAMGRID Plot beam gain surfaces in directed-cosine space.
%
% Usage:
%   plot_beamgrid(BeamGrid)

    if nargin < 1 || isempty(BeamGrid)
        error('plot_beamgrid requires a BeamGrid struct.');
    end

    out_dir = fullfile(pwd, 'figures');
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    beam_subdir = fullfile(out_dir, 'beamgrid');
    if ~isfolder(beam_subdir)
        mkdir(beam_subdir);
    end

    num_samples = 121;
    visible_state = 'off';
    [u_axis, v_axis] = build_uv_axes(BeamGrid, num_samples);
    [U, V] = meshgrid(u_axis, v_axis);
    inside_unit_disk = (U.^2 + V.^2) <= 1;

    gain_fn = get_gain_function(BeamGrid);
    num_beams = get_beam_count(BeamGrid);

    figure_handles = struct();

    overview_fig = figure('Color', 'w', 'Name', 'BeamGrid Overview', 'Visible', visible_state);
    ax = axes(overview_fig);

    layout_name = safe_name(get_beamgrid_label(BeamGrid));
    
    max_gain_linear = zeros(size(U));
    for beam_idx = 1:num_beams
        beam_center_u = BeamGrid.b_u(beam_idx);
        beam_center_v = BeamGrid.b_v(beam_idx);

        gain_linear = gain_fn(U - beam_center_u, V - beam_center_v, beam_idx);
        max_gain_linear = max(max_gain_linear, gain_linear);
    end

    gain_dB = 10 * log10(max(max_gain_linear, eps));
    gain_dB(~inside_unit_disk) = NaN;

    surf(ax, U, V, gain_dB, 'EdgeColor', 'none');
    hold(ax, 'on');
    max_val = max(gain_dB(:), [], 'omitnan');
    plot3(ax, BeamGrid.b_u, BeamGrid.b_v, max_val * ones(size(BeamGrid.b_u)), 'k.', 'MarkerSize', 16);
    hold(ax, 'off');
    view(ax, 45, 35);
    axis(ax, 'tight');
    xlabel(ax, 'u');
    ylabel(ax, 'v');
    zlabel(ax, 'Gain (dB)');
    title(ax, 'Combined Beam Pattern (Max Gain)');
    grid(ax, 'on');
    colormap(ax, turbo);

    exportgraphics(overview_fig, fullfile(beam_subdir, sprintf('%s_beam_surfaces.png', layout_name)), 'Resolution', 300);
    figure_handles.overview = overview_fig;
end

function gain_fn = get_gain_function(BeamGrid)
    if isfield(BeamGrid, 'BeamModel') && isstruct(BeamGrid.BeamModel) && isfield(BeamGrid.BeamModel, 'GainAtOffset') && ~isempty(BeamGrid.BeamModel.GainAtOffset)
        gain_fn = @(du, dv, ~) BeamGrid.BeamModel.GainAtOffset(du, dv);
        return;
    end

    r_u = get_radius(BeamGrid, 'r_beam_u');
    r_v = get_radius(BeamGrid, 'r_beam_v');
    r = get_radius(BeamGrid, 'r_beam');

    gain_fn = @(du, dv, ~) fallback_gaussian_gain(du, dv, r, r_u, r_v);
end

function gain_linear = fallback_gaussian_gain(du, dv, r, r_u, r_v)
    if r_u ~= r_v
        gain_linear = exp(-log(2) * ((du ./ max(r_u, eps)).^2 + (dv ./ max(r_v, eps)).^2));
    else
        gain_linear = exp(-log(2) * ((du.^2 + dv.^2) ./ max(r, eps).^2));
    end
end

function radius = get_radius(BeamGrid, field_name)
    if isfield(BeamGrid, field_name) && ~isempty(BeamGrid.(field_name))
        radius = BeamGrid.(field_name);
    else
        radius = 1;
    end
end

function [u_axis, v_axis] = build_uv_axes(BeamGrid, num_samples)
    if isfield(BeamGrid, 'b_u') && ~isempty(BeamGrid.b_u)
        b_u = BeamGrid.b_u(:);
    else
        b_u = 0;
    end

    if isfield(BeamGrid, 'b_v') && ~isempty(BeamGrid.b_v)
        b_v = BeamGrid.b_v(:);
    else
        b_v = 0;
    end

    r_u = get_radius(BeamGrid, 'r_beam_u');
    r_v = get_radius(BeamGrid, 'r_beam_v');

    u_min = max(-1, min(b_u) - 1.5 * r_u);
    u_max = min(1, max(b_u) + 1.5 * r_u);
    v_min = max(-1, min(b_v) - 1.5 * r_v);
    v_max = min(1, max(b_v) + 1.5 * r_v);

    if ~(isfinite(u_min) && isfinite(u_max) && u_max > u_min)
        u_min = -1;
        u_max = 1;
    end
    if ~(isfinite(v_min) && isfinite(v_max) && v_max > v_min)
        v_min = -1;
        v_max = 1;
    end

    u_axis = linspace(u_min, u_max, num_samples);
    v_axis = linspace(v_min, v_max, num_samples);
end

function num_beams = get_beam_count(BeamGrid)
    if isfield(BeamGrid, 'num_beams') && ~isempty(BeamGrid.num_beams)
        num_beams = BeamGrid.num_beams;
    elseif isfield(BeamGrid, 'b_u')
        num_beams = numel(BeamGrid.b_u);
    else
        num_beams = 0;
    end
end

function nrows = best_tiling(num_beams)
    nrows = max(1, floor(sqrt(max(num_beams, 1))));
end

function ncols = best_tiling_columns(num_beams)
    ncols = max(1, ceil(max(num_beams, 1) / best_tiling(num_beams)));
end

function label = get_beamgrid_label(BeamGrid)
    if isfield(BeamGrid, 'layout') && ~isempty(BeamGrid.layout)
        label = char(BeamGrid.layout);
    else
        label = 'beamgrid';
    end
end

function safe = safe_name(label)
    safe = regexprep(lower(string(label)), '[^a-z0-9_-]+', '_');
    safe = char(strtrim(safe));
    if isempty(safe)
        safe = 'beamgrid';
    end
end