function o = defaultFusionSettings()
    o = struct();
    o.neighborhood_mode = 8;
    o.unreliable_score_strength = 0.0;
    o.unreliable_score = 0;
    o.min_potential = 0.001;
    o.threshold_saturation_gap = 0.05;
    o.discard_unreliable_maps = 1;
    o.cand_map_filter_threshold = 0.1;
    o.verbose = false;
end