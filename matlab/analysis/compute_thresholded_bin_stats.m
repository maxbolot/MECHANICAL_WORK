function stats = compute_thresholded_bin_stats(cfg)
% COMPUTE_THRESHOLDED_BIN_STATS  Compute per-percentile-rank-bin averages of
%   work and lift for a given scenario/region configuration.
%
%   stats = compute_thresholded_bin_stats(cfg)
%
%   cfg must contain:
%     cfg.scenario  — from scenario_control() / scenario_plus4k()
%     cfg.region    — from region_tropics() / etc.
%
%   Returned struct fields:
%     bin_lower, bin_upper, bin_mass, bin_center  (nbin x 1)
%     work_bin_avg, lift_bin_avg, ratio_bin        (nbin x 1)
%     labels                                        (nbin x 1 cell)

    ncfile = cfg.scenario.thresholded_ncfile;
    ensure_input_exists(ncfile);

    [lat, lon, time, percentile] = read_thresholded_coordinates(ncfile);

    mask = build_region_mask(lat, lon, cfg.region);

    work4d = read_thresholded_var(ncfile, 'work');
    lift4d = read_thresholded_var(ncfile, 'lift');

    ntime = size(work4d, 3);
    nperc = size(work4d, 4);

    work_region = work4d(mask.lon_idx, mask.lat_idx, :, :);
    lift_region = lift4d(mask.lon_idx, mask.lat_idx, :, :);
    weight_4d = reshape(mask.lat_weights, 1, [], 1, 1);

    work_num = squeeze(sum(work_region .* weight_4d, [1, 2], 'omitnan'));
    work_den = squeeze(sum(double(~isnan(work_region)) .* weight_4d, [1, 2], 'omitnan'));
    work_spatial_avg = orient_time_percentile(work_num ./ work_den, ntime, nperc);

    lift_num = squeeze(sum(lift_region .* weight_4d, [1, 2], 'omitnan'));
    lift_den = squeeze(sum(double(~isnan(lift_region)) .* weight_4d, [1, 2], 'omitnan'));
    lift_spatial_avg = orient_time_percentile(lift_num ./ lift_den, ntime, nperc);

    [time_weights_days, ~] = cfg.scenario.time_weight_fn(time, ncfile);

    work_avg = nan(nperc, 1);
    lift_avg = nan(nperc, 1);
    for ip = 1:nperc
        work_avg(ip) = weighted_nanmean(work_spatial_avg(:, ip), time_weights_days);
        lift_avg(ip) = weighted_nanmean(lift_spatial_avg(:, ip), time_weights_days);
    end

    [work_total_avg, lift_total_avg] = compute_unthresholded_region_time_means(cfg, mask);

    stats = build_percentile_bin_stats(percentile, work_avg, lift_avg, work_total_avg, lift_total_avg);
end

% -------------------------------------------------------------------------
% Local helpers (mirrors of those in run_thresholded_region_analysis.m)
% -------------------------------------------------------------------------

function ensure_input_exists(ncfile)
    if ~isfile(ncfile)
        error('File not found: %s', ncfile);
    end
end

function [lat, lon, time, percentile] = read_thresholded_coordinates(ncfile)
    try
        lat        = ncread(ncfile, 'lat');
        lon        = ncread(ncfile, 'lon');
        time       = ncread(ncfile, 'time');
        percentile = ncread(ncfile, 'percentile');
    catch ME
        error('Error reading dimensions: %s', ME.message);
    end
end

function values = read_thresholded_var(ncfile, varname)
    try
        raw_values = ncread(ncfile, varname);
        dim_names  = get_var_dim_names(ncfile, varname);
        values     = permute_to_lon_lat_time_percentile(raw_values, dim_names);
    catch ME
        error('Error reading %s variable: %s', varname, ME.message);
    end
end

function dim_names = get_var_dim_names(ncfile, varname)
    info = ncinfo(ncfile, varname);
    dim_names = strings(numel(info.Dimensions), 1);
    for i = 1:numel(info.Dimensions)
        dim_names(i) = string(info.Dimensions(i).Name);
    end
end

function values = permute_to_lon_lat_time_percentile(values_in, dim_names)
    desired    = ["lon", "lat", "time", "percentile"];
    perm_order = zeros(1, numel(desired));
    for i = 1:numel(desired)
        idx = find(dim_names == desired(i), 1);
        if isempty(idx)
            error('Variable is missing required dimension "%s".', desired(i));
        end
        perm_order(i) = idx;
    end
    values = permute(values_in, perm_order);
end

function values = orient_time_percentile(values_in, ntime, nperc)
    values = values_in;
    if size(values, 1) ~= ntime && size(values, 2) == ntime
        values = values.';
    end
    if size(values, 1) ~= ntime || size(values, 2) ~= nperc
        error('Unexpected averaged field shape: got=%s expected=[%d %d]', mat2str(size(values)), ntime, nperc);
    end
end

function [work_avg, lift_avg] = compute_unthresholded_region_time_means(cfg, mask)
    ncfile = cfg.scenario.standard_ncfile;
    ensure_input_exists(ncfile);

    try
        lat_s  = ncread(ncfile, 'lat');  %#ok<NASGU>
        lon_s  = ncread(ncfile, 'lon');  %#ok<NASGU>
        time_s = ncread(ncfile, 'time');
    catch ME
        error('Error reading standard coordinates: %s', ME.message);
    end

    try
        work3d = ncread(ncfile, 'work');
        lift3d = ncread(ncfile, 'lift');
    catch ME
        error('Error reading standard variables: %s', ME.message);
    end

    work_region = work3d(mask.lon_idx, mask.lat_idx, :);
    lift_region = lift3d(mask.lon_idx, mask.lat_idx, :);
    weight_3d   = reshape(mask.lat_weights, 1, [], 1);

    work_num = sum(work_region .* weight_3d, [1, 2], 'omitnan');
    work_den = sum(double(~isnan(work_region)) .* weight_3d, [1, 2], 'omitnan');
    work_spatial_avg = squeeze(work_num ./ work_den);

    lift_num = sum(lift_region .* weight_3d, [1, 2], 'omitnan');
    lift_den = sum(double(~isnan(lift_region)) .* weight_3d, [1, 2], 'omitnan');
    lift_spatial_avg = squeeze(lift_num ./ lift_den);

    time_weights_days = cfg.scenario.time_weight_fn(time_s, ncfile);
    if iscell(time_weights_days)
        time_weights_days = time_weights_days{1};
    end

    work_avg = weighted_nanmean(work_spatial_avg, time_weights_days);
    lift_avg = weighted_nanmean(lift_spatial_avg, time_weights_days);
end

function stats = build_percentile_bin_stats(percentile, work_tail_avg, lift_tail_avg, work_total_avg, lift_total_avg)
    p     = double(percentile(:));
    nperc = numel(p);

    lower = [0.0; p];
    upper = [p; 1.0];
    mass  = upper - lower;
    nbin  = numel(mass);

    work_bin_avg = nan(nbin, 1);
    lift_bin_avg = nan(nbin, 1);

    tail_mass_first = 1.0 - p(1);
    work_bin_avg(1) = (work_total_avg - tail_mass_first * work_tail_avg(1)) / mass(1);
    lift_bin_avg(1) = (lift_total_avg - tail_mass_first * lift_tail_avg(1)) / mass(1);

    for ib = 2:nperc
        p_lo         = p(ib - 1);
        p_hi         = p(ib);
        tail_mass_lo = 1.0 - p_lo;
        tail_mass_hi = 1.0 - p_hi;
        work_bin_avg(ib) = (tail_mass_lo * work_tail_avg(ib - 1) - tail_mass_hi * work_tail_avg(ib)) / mass(ib);
        lift_bin_avg(ib) = (tail_mass_lo * lift_tail_avg(ib - 1) - tail_mass_hi * lift_tail_avg(ib)) / mass(ib);
    end

    work_bin_avg(end) = work_tail_avg(end);
    lift_bin_avg(end) = lift_tail_avg(end);

    ratio_bin = lift_bin_avg ./ work_bin_avg;
    ratio_bin(abs(work_bin_avg) < eps) = NaN;

    labels    = cell(nbin, 1);
    labels{1} = sprintf('<P%s', format_percentile_value(100.0 * p(1)));
    for ib = 2:nperc
        labels{ib} = sprintf('P%s-P%s', format_percentile_value(100.0 * p(ib - 1)), format_percentile_value(100.0 * p(ib)));
    end
    labels{end} = sprintf('>P%s', format_percentile_value(100.0 * p(end)));

    stats            = struct();
    stats.bin_lower  = lower;
    stats.bin_upper  = upper;
    stats.bin_mass   = mass;
    stats.bin_center = 0.5 * (lower + upper);
    stats.work_bin_avg = work_bin_avg;
    stats.lift_bin_avg = lift_bin_avg;
    stats.ratio_bin  = ratio_bin;
    stats.labels     = labels;
end

function text_out = format_percentile_value(value_in)
    text_out = regexprep(sprintf('%.6f', value_in), '0+$', '');
    text_out = regexprep(text_out, '\.$', '');
end
