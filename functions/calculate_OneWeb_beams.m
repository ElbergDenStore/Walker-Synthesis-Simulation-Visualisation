function BeamGrid = calculate_OneWeb_beams(f, G_tx, orbit_height, Min_Elev_deg, FRF)
% CALCULATE_ONEWEB_BEAMS
% Create a OneWeb-style 16-beam fan-beam layout as a custom beam grid.
% The default generic beam grid remains elsewhere; this helper is opt-in.
%
% Model assumption:
% - 16 beams arranged in a 16x1 row
% - each beam is elliptical due to a 32x1 phased-array aperture
% - the 3 dB contours touch along the narrow axis

    if nargin < 4, Min_Elev_deg = 20; end
    if nargin < 5, FRF = 1; end

    Re = 6371; % Earth radius in km
    h = orbit_height / 1000;
    c = physconst('LightSpeed');
    lambda = c / f;

    % Start from the gain-based equivalent beamwidth, then stretch it into an
    % ellipse using the 32x1 element aspect ratio.
    Beamwidth_eq_deg = sqrt(32400 ./ (10.^(G_tx/10)));
    element_dims = [32, 1];
    aspect_ratio = element_dims(1) / element_dims(2);

    % Keep the same equivalent directivity, but shape it as an ellipse.
    Beamwidth_minor_deg = Beamwidth_eq_deg / sqrt(aspect_ratio);
    Beamwidth_major_deg = Beamwidth_eq_deg * sqrt(aspect_ratio);

    r_beam_minor = sind(Beamwidth_minor_deg / 2);
    r_beam_major = sind(min(Beamwidth_major_deg / 2, 89.9));
    eta_max = asind((Re / (Re + h)) * cosd(Min_Elev_deg));

    % 16 beams in a single row. The narrow axis is along u, so the 3 dB
    % contours just touch from beam to beam along the row.
    num_beams = 16;
    center_spacing_u = 2 * r_beam_minor;
    row_offsets = (0:num_beams-1) - (num_beams - 1) / 2;
    b_u = row_offsets(:) * center_spacing_u;
    b_v = zeros(num_beams, 1);

    % Metadata for compatibility with the existing codebase.
    b_q = (0:num_beams-1)';
    b_r = zeros(num_beams, 1);

    % Neighbor indices: immediate left/right are the only direct neighbors in
    % this fan-beam arrangement. Remaining slots stay empty.
    neighbor_idx = nan(num_beams, 6);
    for b = 1:num_beams
        if b > 1
            neighbor_idx(b, 1) = b - 1;
        end
        if b < num_beams
            neighbor_idx(b, 2) = b + 1;
        end
    end

    % Frequency reuse group assignment: 8-channel pattern repeated over 16 beams.
    channel_id = mod((0:num_beams-1)', 8) + 1;

    BeamGrid.b_u = b_u;
    BeamGrid.b_v = b_v;
    BeamGrid.b_q = b_q;
    BeamGrid.b_r = b_r;
    BeamGrid.neighbor_idx = neighbor_idx;
    BeamGrid.r_beam = r_beam_minor;
    BeamGrid.r_beam_u = r_beam_minor;
    BeamGrid.r_beam_v = r_beam_major;
    BeamGrid.r_beam_minor = r_beam_minor;
    BeamGrid.r_beam_major = r_beam_major;
    BeamGrid.lambda = lambda;
    BeamGrid.Beamwidth_deg = Beamwidth_eq_deg;
    BeamGrid.Beamwidth_minor_deg = Beamwidth_minor_deg;
    BeamGrid.Beamwidth_major_deg = Beamwidth_major_deg;
    BeamGrid.num_beams = num_beams;
    BeamGrid.FRF = FRF;
    BeamGrid.channel_id = channel_id;
    BeamGrid.layout = 'OneWeb-16x1-elliptical';
    BeamGrid.element_dims = element_dims;
    BeamGrid.beam_shape = 'elliptical';
    BeamGrid.eta_max = eta_max;
end
