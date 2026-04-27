function out = compute_control_work_lift_latband_bar_chart_stats(cfg)
    % Computes control-only work/lift regional summaries for mean and thresholded conditions.
    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    if ~isfield(cfg, 'scenario') || isempty(cfg.scenario)
        cfg.scenario = scenario_control();
    end
    if ~isfield(cfg, 'threshold_percentile') || ~isfinite(cfg.threshold_percentile)
        cfg.threshold_percentile = 0.99999;
    end
    if ~isfield(cfg, 'regions') || isempty(cfg.regions)
        cfg.regions = default_regions_local();
    end

    ensure_input_exists(cfg.scenario.standard_ncfile);
    ensure_input_exists(cfg.scenario.thresholded_by_lat_band_ncfile);

    [lat_std, lon_std, time_std] = read_standard_coordinates(cfg.scenario.standard_ncfile);
    work_std = read_required_var(cfg.scenario.standard_ncfile, 'work');
    lift_std = read_required_var(cfg.scenario.standard_ncfile, 'lift');
    [time_weights_std, missing_steps_std] = cfg.scenario.time_weight_fn(time_std, cfg.scenario.standard_ncfile);

    [lat_thr, lon_thr, time_thr] = read_standard_coordinates(cfg.scenario.thresholded_by_lat_band_ncfile);
    [work_thr, selected_percentile] = read_lon_lat_time_from_var(cfg.scenario.thresholded_by_lat_band_ncfile, 'work', cfg.threshold_percentile);
    [lift_thr, ~] = read_lon_lat_time_from_var(cfg.scenario.thresholded_by_lat_band_ncfile, 'lift', cfg.threshold_percentile);
    [event_count, ~] = read_lon_lat_time_from_var(cfg.scenario.thresholded_by_lat_band_ncfile, 'event_count', cfg.threshold_percentile);
    [time_weights_thr, missing_steps_thr] = cfg.scenario.time_weight_fn(time_thr, cfg.scenario.thresholded_by_lat_band_ncfile);

    validate_compatible_grid_local(lat_std, lon_std, lat_thr, lon_thr);

    nreg = numel(cfg.regions);
    rows = repmat(struct( ...
        'name', '', ...
        'mean_work', NaN, ...
        'mean_lift', NaN, ...
        'threshold_work', NaN, ...
        'threshold_lift', NaN, ...
        'mean_ratio', NaN, ...
        'threshold_ratio', NaN, ...
        'mean_missing_steps', missing_steps_std, ...
        'threshold_missing_steps', missing_steps_thr, ...
        'selected_percentile', selected_percentile), nreg, 1);

    for ireg = 1:nreg
        region = cfg.regions(ireg);
        mask_std = build_latband_mask_local(lat_std, lon_std, region);
        mask_thr = build_latband_mask_local(lat_thr, lon_thr, region);

        mean_work = weighted_region_mean_local(work_std, mask_std, time_weights_std, []);
        mean_lift = weighted_region_mean_local(lift_std, mask_std, time_weights_std, []);
        threshold_work = weighted_region_mean_local(work_thr, mask_thr, time_weights_thr, event_count);
        threshold_lift = weighted_region_mean_local(lift_thr, mask_thr, time_weights_thr, event_count);

        rows(ireg).name = region.name;
        rows(ireg).mean_work = mean_work;
        rows(ireg).mean_lift = mean_lift;
        rows(ireg).threshold_work = threshold_work;
        rows(ireg).threshold_lift = threshold_lift;
        rows(ireg).mean_ratio = safe_ratio_local(mean_lift, mean_work);
        rows(ireg).threshold_ratio = safe_ratio_local(threshold_lift, threshold_work);
    end

    out = struct();
    out.scenario = cfg.scenario;
    out.threshold_percentile_requested = cfg.threshold_percentile;
    out.threshold_percentile_selected = selected_percentile;
    out.regions = cfg.regions;
    out.rows = rows;
end


function regions = default_regions_local()
    regions = [ ...
        struct('name', 'tropics', 'label', 'Tropics', 'lat_ranges', [-30, 30]), ...
        struct('name', 'midlatitudes', 'label', 'Midlatitudes', 'lat_ranges', [-60, -30; 30, 60]) ...
    ];
end


function validate_compatible_grid_local(lat_a, lon_a, lat_b, lon_b)
    if numel(lat_a) ~= numel(lat_b) || any(abs(double(lat_a(:)) - double(lat_b(:))) > 1.0e-10)
        error('Standard and thresholded latitude coordinates do not match.');
    end
    if numel(lon_a) ~= numel(lon_b) || any(abs(double(lon_a(:)) - double(lon_b(:))) > 1.0e-10)
        error('Standard and thresholded longitude coordinates do not match.');
    end
end


function mask = build_latband_mask_local(lat, lon, region)
    lat = double(lat(:));
    lon = double(lon(:));

    lat_mask = false(size(lat));
    for irange = 1:size(region.lat_ranges, 1)
        lat_mask = lat_mask | (lat >= region.lat_ranges(irange, 1) & lat <= region.lat_ranges(irange, 2));
    end

    lat_idx = find(lat_mask);
    lon_idx = 1:numel(lon);
    if isempty(lat_idx)
        error('No latitude points found for region %s.', region.name);
    end

    lat_weights = cosd(lat(lat_idx));
    if all(lat_weights <= 0)
        error('Latitude weights are non-positive for region %s.', region.name);
    end

    mask = struct();
    mask.lat_idx = lat_idx;
    mask.lon_idx = lon_idx;
    mask.lat_weights = lat_weights(:);
end


function avg = weighted_region_mean_local(field3d, mask, time_weights_days, extra_cell_weights)
    subset = double(field3d(mask.lon_idx, mask.lat_idx, :));
    lat_weights = reshape(mask.lat_weights, 1, [], 1);
    time_weights = reshape(double(time_weights_days(:)), 1, 1, []);
    weights3d = repmat(lat_weights .* time_weights, size(subset, 1), 1, 1);

    if nargin >= 4 && ~isempty(extra_cell_weights)
        extra_subset = double(extra_cell_weights(mask.lon_idx, mask.lat_idx, :));
        weights3d = weights3d .* extra_subset;
    end

    valid = isfinite(subset) & isfinite(weights3d) & (weights3d > 0);
    if ~any(valid(:))
        avg = NaN;
        return;
    end

    avg = sum(subset(valid) .* weights3d(valid)) / sum(weights3d(valid));
end


function ratio = safe_ratio_local(numerator, denominator)
    if ~isfinite(numerator) || ~isfinite(denominator) || abs(denominator) < eps
        ratio = NaN;
        return;
    end
    ratio = numerator / denominator;
end