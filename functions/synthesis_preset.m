function [Master_config, default_heights_km] = synthesis_preset(preset)
% SYNTHESIS_PRESET  Named configuration blocks for numerical_walker_synthesis.
%
%   [Master_config, default_heights_km] = synthesis_preset(preset)
%
%   Returns the default Master_config and altitude vector for a named
%   synthesis scenario. The calling script (numerical_walker_synthesis.m)
%   may override individual fields after this call.
%
%   preset (string):
%     "regional_delta" - Walker Delta, Arctic/Nordic region (54.6-83.7 N),
%                        inclination swept 70-80 deg.
%     "global_delta"   - Walker Delta, equatorial/tropical band (0-30 N),
%                        wider plane/sat search.
%     "star"           - Walker Star, polar (i = 90 deg), Arctic/Nordic region.
%
%   The statistical "Detailed" tier (Duration/Num_UEs) is derived from a
%   coverage-confidence sample-count formula, identical to the original
%   Find_valid_* scripts.

    arguments
        preset (1,1) string
    end

    cfg = struct();

    switch lower(preset)
        case "regional_delta"
            default_heights_km          = 500:10:1200;
            cfg.WalkerStar              = false;
            cfg.folder_prefix           = "Master_Sweep";
            cfg.Lat_range_deg           = [54+(35/60), 83+(40/60)];  % Denmark min, Greenland max
            cfg.Min_elevation_UE        = 20;
            cfg.Num_Planes              = 2:25;
            cfg.Sats_Plane              = 2:25;
            cfg.Inc_vec                 = linspace(70, 80, 111);
            cfg.Target_num_candidates   = 1;
            cfg.SampleTime              = 660;
            cfg.Ultrafast.Duration_h    = 3;
            cfg.Ultrafast.Num_UEs       = 200;
            cfg.Fast.Duration_h         = 50;
            cfg.Fast.Num_UEs            = 800;
            cfg.Max_workers             = 16;
            certainty = 0.99; fractional_area = 0.001; fractional_time = 0.001;

        case "iris2"
            default_heights_km          = [1200 8000];
            cfg.WalkerStar              = false;
            cfg.folder_prefix           = "Iris2";
            cfg.Lat_range_deg           = [30, 60];  % Denmark min, Greenland max
            cfg.Min_elevation_UE        = 20;
            cfg.Num_Planes              = 2:25;
            cfg.Sats_Plane              = 2:40;
            cfg.Inc_vec                 = linspace(20, 80, 61);
            cfg.Target_num_candidates   = 20;
            cfg.SampleTime              = 660;
            cfg.Ultrafast.Duration_h    = 3;
            cfg.Ultrafast.Num_UEs       = 200;
            cfg.Fast.Duration_h         = 50;
            cfg.Fast.Num_UEs            = 800;
            cfg.Max_workers             = 16;
            certainty = 0.99; fractional_area = 0.001; fractional_time = 0.001;

        case "global_delta"
            default_heights_km          = 500:100:1200;
            cfg.WalkerStar              = false;
            cfg.Lat_range_deg           = [0, 30];
            cfg.folder_prefix           = sprintf("Master_Sweep_lat%.0f_%.0f", ...
                                                  cfg.Lat_range_deg(1), cfg.Lat_range_deg(2));
            cfg.Min_elevation_UE        = 20;
            cfg.Num_Planes              = 1:25;
            cfg.Sats_Plane              = 4:40;
            cfg.Inc_vec                 = linspace(max(cfg.Lat_range_deg)-15, ...
                                                   min(max(cfg.Lat_range_deg), 80), 16);
            cfg.Target_num_candidates   = 1;
            cfg.SampleTime              = 660;
            cfg.Ultrafast.Duration_h    = 1;
            cfg.Ultrafast.Num_UEs       = 1000;
            cfg.Fast.Duration_h         = 50;
            cfg.Fast.Num_UEs            = 2000;
            cfg.Max_workers             = 16;
            certainty = 0.99; fractional_area = 0.0005; fractional_time = 0.001;

        case "star"
            default_heights_km          = 500:10:1200;
            cfg.WalkerStar              = true;
            cfg.folder_prefix           = "Master_Sweep_WalkerStar";
            cfg.Lat_range_deg           = [54+(35/60), 83+(40/60)];
            cfg.Min_elevation_UE        = 20;
            cfg.Num_Planes              = 2:15;
            cfg.Sats_Plane              = 5:35;
            cfg.Inc_vec                 = 90;
            cfg.Target_num_candidates   = 1;
            cfg.SampleTime              = 660;
            cfg.Ultrafast.Duration_h    = 3;
            cfg.Ultrafast.Num_UEs       = 200;
            cfg.Fast.Duration_h         = 50;
            cfg.Fast.Num_UEs            = 800;
            cfg.Worker_stall_timeout_s  = 3600;   % >> longest legitimate detailed run
            cfg.Max_workers             = 4;
            certainty = 0.99; fractional_area = 0.001; fractional_time = 0.001;

        otherwise
            error("synthesis_preset:UnknownPreset", ...
                "Unknown preset '%s'. Use 'regional_delta', 'global_delta', or 'star'.", preset);
    end

    % --- Derived "Detailed" statistical tier (coverage-confidence sampling) ---
    required_samples       = log(1 - certainty) / log(1 - fractional_area * fractional_time);
    cfg.Detailed.Duration_h = ceil(sqrt(required_samples)) / (3600 / cfg.SampleTime);
    cfg.Detailed.Num_UEs    = ceil(sqrt(required_samples));

    Master_config = cfg;
end
