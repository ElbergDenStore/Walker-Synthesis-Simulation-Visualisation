% PLOT_SEAM_GEOMETRY  Napier's Circle derivation and counter-rotating seam geometry.
%   Produces five PNG figures:
%     napiers_circle.png           – annotated pentagon showing the SIN-TAAD derivation
%     seam_longitude_drift.png     – per-satellite longitude drift Δλ(φ)
%     seam_gap_vs_latitude.png     – seam longitude gap growing with latitude
%     seam_track_bowing.png        – counter-rotating track bowing visualisation
%     seam_gc_gap_vs_latitude.png  – great-circle seam gap vs latitude
%
% Counter-rotating seam geometry summary:
%   Ascending track (plane 1, RAAN=0): at latitude φ drifts EAST by Δλ(φ)
%   Descending track (plane P): at latitude φ drifts WEST by Δλ(φ) from its desc. node
%   → Seam longitude gap = RAAN_seam + 2·Δλ(φ)   (grows with latitude)
%   → Physical gap = arccos(sin²φ + cos²φ·cos(lon_gap))  (peaks at intermediate lat)
%
% Usage:
%   plot_seam_geometry()                    – saves to plotting_scripts/figures/
%   plot_seam_geometry('path/to/output')    – saves to the given directory

save_path = '';   % output folder for the seam-geometry PNGs ('' = figures/seam_gap)

script_dir = fileparts(mfilename('fullpath'));
if isempty(save_path)
    save_path = fullfile(script_dir, '../figures/seam_gap');
end

set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

%% ── Figure 1: Napier's Circle ────────────────────────────────────────────
f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 400], ...
    'Name', 'Napier''s Circle');
ax1 = axes('Parent', f1);
hold(ax1, 'on');
axis(ax1, 'equal');
axis(ax1, 'off');

% Five sectors of Napier's Pentagon (72 degrees each), starting from top
angles = linspace(pi/2, pi/2 - 2*pi, 6);

% Pentagon labels: a, b, co-A, co-c, co-B
labels = {
    'Latitude (a) = \phi', ...
    'Drift (b) = \Delta\lambda', ...
    'co-Inclination (co-A) = 90^\circ - i', ...
    'co-Distance (co-c) = 90^\circ - d', ...
    'co-Angle (co-B) = 90^\circ - \beta'
};

% SIN-TAAD roles (highlight the active triplet)
roles  = {'ADJACENT 1', 'MIDDLE PART', 'ADJACENT 2', '', ''};

% Yellow = adjacent, green = middle, grey = unused
colors = {[1 0.85 0.4], [0.4 0.85 0.4], [1 0.85 0.4], [0.9 0.9 0.9], [0.9 0.9 0.9]};

for k = 1:5
    th_fill = linspace(angles(k), angles(k+1), 50);
    fill(ax1, [0 cos(th_fill) 0], [0 sin(th_fill) 0], colors{k}, ...
        'EdgeColor', 'k', 'LineWidth', 2);
    plot(ax1, [0 cos(angles(k))], [0 sin(angles(k))], 'k', 'LineWidth', 2);

    mid_angle = (angles(k) + angles(k+1)) / 2;

    text(0.65*cos(mid_angle), 0.65*sin(mid_angle), labels{k}, ...
        'Parent', ax1, 'HorizontalAlignment', 'center', ...
        'FontSize', 11, 'FontWeight', 'bold');

    if ~isempty(roles{k})
        text(0.35*cos(mid_angle), 0.35*sin(mid_angle), roles{k}, ...
            'Parent', ax1, 'HorizontalAlignment', 'center', ...
            'FontSize', 9, 'Color', 'k', 'FontAngle', 'italic', 'FontWeight', 'bold');
    end
end

title(ax1, 'Napier''s Circle: Solving for Longitude Drift', ...
    'FontSize', 15, 'FontWeight', 'bold');

text(0, -1.20, 'The SIN-TAAD Rule:', ...
    'Parent', ax1, 'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
text(0, -1.36, 'sin(MIDDLE) = tan(ADJACENT 1) \times tan(ADJACENT 2)', ...
    'Parent', ax1, 'HorizontalAlignment', 'center', 'FontSize', 12);
text(0, -1.55, 'sin(\Delta\lambda) = tan(\phi) \times tan(90^\circ - i)', ...
    'Parent', ax1, 'HorizontalAlignment', 'center', ...
    'FontSize', 13, 'Color', '#D95319', 'FontWeight', 'bold');
text(0, -1.73, '\Downarrow', ...
    'Parent', ax1, 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
text(0, -1.90, '\Delta\lambda = arcsin( tan(\phi) / tan(i) )', ...
    'Parent', ax1, 'HorizontalAlignment', 'center', ...
    'FontSize', 13, 'Color', '#0072BD', 'FontWeight', 'bold');

xlim(ax1, [-1.5 1.5]);
ylim(ax1, [-2.1 1.2]);

out1 = fullfile(save_path, 'napiers_circle.png');
exportgraphics(f1, out1, 'Resolution', 300);
close(f1);
fprintf('Saved: %s\n', out1);

%% ── Figure 2a: Per-satellite longitude drift ────────────────────────────
% From Napier's circle:  Δλ(φ) = arcsin( tan(φ) / tan(i) )
%
% Ascending track  (plane 1, RAAN = 0):          lon(φ) = +Δλ(φ)  (bows east)
% Descending track (plane P, desc. node at
%   360°−RAAN_seam):                              lon(φ) = (360°−RAAN_seam) − Δλ(φ)
%                                                          (bows west)
%
% Short-way seam longitude gap:
%   gap(φ) = 360° − [(360°−RAAN_seam) − Δλ] + Δλ
%          = RAAN_seam + 2·Δλ(φ)   → gap GROWS from RAAN_seam to RAAN_seam+180°

inc_array    = [80, 87, 90];
seam_gap_deg = 26;   % equatorial RAAN gap between the two seam planes (degrees)
line_colors  = {'#0072BD', '#D95319', '#EDB120'};

f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 400], ...
    'Name', 'Seam Gap Geometry');

%% ── Left panel: Δλ(φ) — per-satellite longitude drift from ascending node
ax2L = subplot(1, 2, 1, 'Parent', f2);
hold(ax2L, 'on');

for i = 1:length(inc_array)
    inc = inc_array(i);
    lat = 0:0.1:inc;
    d_lon = asind(tand(lat) ./ tand(inc));   % Napier's circle: Δλ(φ)
    plot(ax2L, lat, d_lon, 'LineWidth', 2, 'Color', line_colors{i}, ...
        'DisplayName', sprintf('%d^\\circ', inc_array(i)));
end

% Reference line: Δλ = 90° is the physical maximum (satellite at orbit apex)
yline(ax2L, 90, 'k--', 'LineWidth', 1.5, 'DisplayName', '\Delta\lambda = 90° (orbit apex)');

legend(ax2L, 'Location', 'northwest', 'FontSize', 11);
xlabel(ax2L, 'Latitude \phi (deg)', 'FontWeight', 'bold');
ylabel(ax2L, '\Delta\lambda (deg)', 'FontWeight', 'bold');
title(ax2L, {'Per-Satellite Longitude Drift', '\Delta\lambda = arcsin( tan(\phi) / tan(i) )'}, 'FontSize', 12);
xlim(ax2L, [0 90]);
grid(ax2L, 'on');
set(ax2L, 'FontSize', 12);

%% ── Right panel: seam longitude gap growing with latitude
ax2R = subplot(1, 2, 2, 'Parent', f2);
hold(ax2R, 'on');

for i = 1:length(inc_array)
    inc = inc_array(i);
    lat = 0:0.01:inc;
    d_lon    = asind(tand(lat) ./ tand(inc));     % per-satellite longitude drift
    lon_gap  = seam_gap_deg + 2 .* d_lon;         % gap GROWS: both tracks bow outward

    plot(ax2R, lat, lon_gap, 'LineWidth', 2, 'Color', line_colors{i}, ...
        'DisplayName', sprintf('%d^\\circ', inc_array(i)));
end

yline(ax2R, 180, 'k--', 'LineWidth', 1.2, 'DisplayName', '180° (antipodal)');
legend(ax2R, 'Location', 'northwest', 'FontSize', 11);
xlabel(ax2R, 'Latitude \phi (deg)', 'FontWeight', 'bold');
ylabel(ax2R, 'Seam Longitude Gap (deg)', 'FontWeight', 'bold');
title(ax2R, {'Seam Longitude Gap Grows with Latitude', ...
    'gap(\phi) = RAAN_{seam} + 2·\Delta\lambda(\phi)'}, 'FontSize', 12);
xlim(ax2R, [0 90]);
ylim(ax2R, [0 seam_gap_deg + 185]);
grid(ax2R, 'on');
set(ax2R, 'FontSize', 12);

out2 = fullfile(save_path, 'seam_gap_vs_latitude.png');
exportgraphics(f2, out2, 'Resolution', 300);
close(f2);
fprintf('Saved: %s\n', out2);

%% ── Figure 3a: Counter-rotating seam track bowing ───────────────────────
% Left:  The two seam tracks drawn as longitude-vs-latitude (centred on seam
%        midpoint).  Both tracks bow OUTWARD as latitude grows — the ascending
%        track drifts east (+Δλ) and the descending track drifts west (−Δλ).
%
% Right: Great-circle arc distance between the two seam tracks.
%        gc_gap(φ) = arccos( sin²φ + cos²φ · cos(RAAN_seam + 2·Δλ(φ)) )
%
%        For near-polar orbits (i→90°), Δλ stays small so gc_gap decreases
%        monotonically — the cosine compression of latitude dominates.
%        For lower inclinations, Δλ grows faster and gc_gap PEAKS at an
%        intermediate latitude before shrinking toward zero near the pole.
%        That peak is the latitude hardest to cover across the seam.

f3a = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 450], ...
    'Name', 'Seam Track Bowing');
ax3a = axes('Parent', f3a);
hold(ax3a, 'on');

for ii = 1:length(inc_array)
    inc  = inc_array(ii);
    lat  = 0:0.1:inc;
    d_lon    = asind(tand(lat) ./ tand(inc));
    asc_lon  =  seam_gap_deg/2 + d_lon;    % ascending track bows east
    desc_lon = -seam_gap_deg/2 - d_lon;    % descending track bows west

    plot(ax3a, asc_lon, lat, '-',  'LineWidth', 2, 'Color', line_colors{ii}, ...
        'DisplayName', sprintf('%d^\\circ', inc_array(ii)));
    plot(ax3a, desc_lon, lat, '--', 'LineWidth', 2, 'Color', line_colors{ii}, ...
        'HandleVisibility', 'off');
end

xline(ax3a, 0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel(ax3a, 'Longitude from seam centre (deg)', 'FontWeight', 'bold');
ylabel(ax3a, 'Latitude (deg)', 'FontWeight', 'bold');
title(ax3a, {'Counter-Rotating Seam: Track Bowing', ...
    'solid = ascending   |   dashed = descending'}, 'FontSize', 12);
legend(ax3a, 'Location', 'northeast', 'FontSize', 11);
ylim(ax3a, [0 90]);
grid(ax3a, 'on');
set(ax3a, 'FontSize', 12);

out3a = fullfile(save_path, 'seam_track_bowing.png');
exportgraphics(f3a, out3a, 'Resolution', 300);
close(f3a);
fprintf('Saved: %s\n', out3a);

%% ── Figure 3b: Great-circle distance between the two seam tracks ─────────
f3b = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 500 450], ...
    'Name', 'Seam Great-Circle Gap');
ax3b = axes('Parent', f3b);
hold(ax3b, 'on');

for ii = 1:length(inc_array)
    inc  = inc_array(ii);
    lat  = 0:0.01:inc;
    d_lon   = asind(tand(lat) ./ tand(inc));
    lon_gap = seam_gap_deg + 2 .* d_lon;   % longitude gap, grows with latitude

    % Spherical law of cosines: arc distance between two points at same latitude
    gc_gap = acosd(sind(lat).^2 + cosd(lat).^2 .* cosd(lon_gap));

    % [gc_max, idx_max] = max(gc_gap);

    plot(ax3b, lat, gc_gap, 'LineWidth', 2, 'Color', line_colors{ii}, ...
        'DisplayName', sprintf('%d^\\circ', ...
        inc_array(ii)));
    % % Mark the peak
    % plot(ax3b, lat(idx_max), gc_max, 'o', 'MarkerSize', 8, ...
    %     'Color', line_colors{ii}, 'MarkerFaceColor', line_colors{ii}, ...
    %     'HandleVisibility', 'off');
end

xlabel(ax3b, 'Latitude (deg)', 'FontWeight', 'bold');
ylabel(ax3b, 'Great-Circle Seam Gap (deg)', 'FontWeight', 'bold');
title(ax3b, 'Physical Distance Between Seam Tracks', 'FontSize', 12);
legend(ax3b, 'Location', 'northeast', 'FontSize', 11);
xlim(ax3b, [0 90]);
grid(ax3b, 'on');
set(ax3b, 'FontSize', 12);

out3b = fullfile(save_path, 'seam_gc_gap_vs_latitude.png');
exportgraphics(f3b, out3b, 'Resolution', 300);
close(f3b);
fprintf('Saved: %s\n', out3b);
