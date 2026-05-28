function out = compute_control_work_lift_latband_bar_chart_panel_a_boxplot_stats(cfg)
    % Computes panel-a means and per-sample distributions for control work/lift by latitude family.
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

    base_stats = compute_control_work_lift_latband_bar_chart_stats(cfg);

    ensure_input_exists(cfg.scenario.standard_ncfile);
    ensure_input_exists(cfg.scenario.thresholded_by_lat_band_ncfile);

    [lat_std, lon_std] = read_standard_coordinates(cfg.scenario.standard_ncfile);
    work_std = read_required_var(cfg.scenario.standard_ncfile, 'work');
    lift_std = read_required_var(cfg.scenario.standard_ncfile, 'lift');

    [lat_thr, lon_thr] = read_standard_coordinates(cfg.scenario.thresholded_by_lat_band_ncfile);
    [work_thr, selected_percentile] = read_lon_lat_time_from_var(cfg.scenario.thresholded_by_lat_band_ncfile, 'work', cfg.threshold_percentile);
    [lift_thr, ~] = read_lon_lat_time_from_var(cfg.scenario.thresholded_by_lat_band_ncfile, 'lift', cfg.threshold_percentile);
    [event_count, ~] = read_lon_lat_time_from_var(cfg.scenario.thresholded_by_lat_band_ncfile, 'event_count', cfg.threshold_percentile);

    validate_compatible_grid_local(lat_std, lon_std, lat_thr, lon_thr);

    rows = base_stats.rows;
    for ireg = 1:numel(cfg.regions)
        region = cfg.regions(ireg);
        mask_std = build_latband_mask_local(lat_std, lon_std, region);
        mask_thr = build_latband_mask_local(lat_thr, lon_thr, region);

        rows(ireg).mean_work_samples = extract_region_samples_local(work_std, mask_std, []);
        rows(ireg).mean_lift_samples = extract_region_samples_local(lift_std, mask_std, []);
        rows(ireg).threshold_work_samples = extract_region_samples_local(work_thr, mask_thr, event_count);
        rows(ireg).threshold_lift_samples = extract_region_samples_local(lift_thr, mask_thr, event_count);
    end

    out = base_stats;
    out.threshold_percentile_selected = selected_percentile;
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

    mask = struct();
    mask.lat_idx = lat_idx;
    mask.lon_idx = lon_idx;
end


function samples = extract_region_samples_local(field3d, mask, event_count)
    subset = double(field3d(mask.lon_idx, mask.lat_idx, :));
    valid = isfinite(subset);

    if nargin >= 3 && ~isempty(event_count)
        evt_subset = double(event_count(mask.lon_idx, mask.lat_idx, :));
        valid = valid & isfinite(evt_subset) & (evt_subset > 0);
    end

    samples = subset(valid);
    samples = samples(:);
end
