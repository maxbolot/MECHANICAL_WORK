function run_thresholded_region_analysis(cfg)
    ncfile = cfg.scenario.thresholded_ncfile;
    ensure_input_exists(ncfile);

    if ~isfield(cfg, 'smooth_before_plotting')
        cfg.smooth_before_plotting = false;
    end
    if ~isfield(cfg, 'smooth_window_days')
        cfg.smooth_window_days = 5;
    end

    fprintf('Reading NetCDF file: %s\n', ncfile);
    [lat, lon, time, percentile] = read_thresholded_coordinates(ncfile);

    mask = build_region_mask(lat, lon, cfg.region);
    fprintf('%s\n', mask.selection_message);

    work4d = read_thresholded_var(ncfile, 'work');
    fprintf('Work variable reordered shape [lon lat time percentile]: %s\n', mat2str(size(work4d)));
    lift4d = read_thresholded_var(ncfile, 'lift');
    fprintf('Lift variable reordered shape [lon lat time percentile]: %s\n', mat2str(size(lift4d)));

    ntime = size(work4d, 3);
    nperc = size(work4d, 4);
    if nperc ~= numel(percentile)
        error('Percentile dimension mismatch: data=%d, coordinate=%d', nperc, numel(percentile));
    end

    work_region = work4d(mask.lon_idx, mask.lat_idx, :, :);
    lift_region = lift4d(mask.lon_idx, mask.lat_idx, :, :);
    weight_4d = reshape(mask.lat_weights, 1, [], 1, 1);

    work_num = squeeze(sum(work_region .* weight_4d, [1, 2], 'omitnan'));
    work_den = squeeze(sum(double(~isnan(work_region)) .* weight_4d, [1, 2], 'omitnan'));
    work_spatial_avg = orient_time_percentile(work_num ./ work_den, ntime, nperc);

    lift_num = squeeze(sum(lift_region .* weight_4d, [1, 2], 'omitnan'));
    lift_den = squeeze(sum(double(~isnan(lift_region)) .* weight_4d, [1, 2], 'omitnan'));
    lift_spatial_avg = orient_time_percentile(lift_num ./ lift_den, ntime, nperc);

    [time_weights_days, missing_steps] = cfg.scenario.time_weight_fn(time, ncfile);
    if numel(time_weights_days) ~= ntime
        error('Time-weight vector length (%d) does not match data time length (%d).', numel(time_weights_days), ntime);
    end

    work_avg = nan(nperc, 1);
    lift_avg = nan(nperc, 1);
    for ip = 1:nperc
        work_avg(ip) = weighted_nanmean(work_spatial_avg(:, ip), time_weights_days);
        lift_avg(ip) = weighted_nanmean(lift_spatial_avg(:, ip), time_weights_days);
    end

    ratio_spatial_avg = lift_spatial_avg ./ work_spatial_avg;
    ratio_spatial_avg(abs(work_spatial_avg) < eps) = NaN;
    ratio_avg = lift_avg ./ work_avg;
    ratio_avg(abs(work_avg) < eps) = NaN;

    [work_total_avg, lift_total_avg] = compute_unthresholded_region_time_means(cfg, mask);
    bin_stats = build_percentile_bin_stats(percentile, work_avg, lift_avg, work_total_avg, lift_total_avg);

    [plot_time, plot_time_label] = build_time_axis(ncfile, time);
    percentile_labels = build_percentile_labels(percentile);

    work_plot = work_spatial_avg;
    lift_plot = lift_spatial_avg;
    ratio_plot = ratio_spatial_avg;
    if cfg.smooth_before_plotting
        if cfg.smooth_window_days < 1 || mod(cfg.smooth_window_days, 1) ~= 0
            error('smooth_window_days must be a positive integer.');
        end
        work_plot = running_mean_omitnan_cols(work_spatial_avg, cfg.smooth_window_days);
        lift_plot = running_mean_omitnan_cols(lift_spatial_avg, cfg.smooth_window_days);
        ratio_plot = lift_plot ./ work_plot;
        ratio_plot(abs(work_plot) < eps) = NaN;
        title_suffix = sprintf(' (All Percentile Thresholds (%d-day running mean))', cfg.smooth_window_days);
    else
        title_suffix = ' (All Percentile Thresholds)';
    end

    summarize_work_lift(struct( ...
        'is_thresholded', true, ...
        'region_label', cfg.region.result_label, ...
        'missing_steps', missing_steps, ...
        'ntime', ntime, ...
        'percentile_labels', {percentile_labels}, ...
        'work_avg', work_avg, ...
        'lift_avg', lift_avg, ...
        'ratio_avg', ratio_avg, ...
        'smooth_before_plotting', cfg.smooth_before_plotting, ...
        'smooth_window_days', cfg.smooth_window_days));

    plot_work_lift_timeseries(struct( ...
        'plot_time', plot_time, ...
        'plot_time_label', plot_time_label, ...
        'work_plot', work_plot, ...
        'lift_plot', lift_plot, ...
        'ratio_plot', ratio_plot, ...
        'title_prefix', cfg.region.title_label, ...
        'title_suffix', title_suffix, ...
        'legend_labels', {percentile_labels}));

    plot_percentile_bin_ratio(bin_stats, cfg.region.title_label);
    print_percentile_bin_summary(bin_stats);
end

function ensure_input_exists(ncfile)
    if ~isfile(ncfile)
        error('File not found: %s', ncfile);
    end
end

function [lat, lon, time, percentile] = read_thresholded_coordinates(ncfile)
    try
        lat = ncread(ncfile, 'lat');
        lon = ncread(ncfile, 'lon');
        time = ncread(ncfile, 'time');
        percentile = ncread(ncfile, 'percentile');
        fprintf('File dimensions: lat=%d, lon=%d, time=%d, percentile=%d\n', length(lat), length(lon), length(time), length(percentile));
    catch ME
        error('Error reading dimensions: %s', ME.message);
    end
end

function values = read_thresholded_var(ncfile, varname)
    try
        raw_values = ncread(ncfile, varname);
        dim_names = get_var_dim_names(ncfile, varname);
        values = permute_to_lon_lat_time_percentile(raw_values, dim_names);
    catch ME
        error('Error reading %s variable: %s\nAvailable variables: %s', varname, ME.message, get_available_vars(ncfile));
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
    desired = ["lon", "lat", "time", "percentile"];
    perm_order = zeros(1, numel(desired));
    for i = 1:numel(desired)
        idx = find(dim_names == desired(i), 1);
        if isempty(idx)
            error('Variable is missing required dimension "%s". Found: %s', desired(i), strjoin(cellstr(dim_names), ', '));
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

function labels = build_percentile_labels(percentile)
    labels = cell(numel(percentile), 1);
    for i = 1:numel(percentile)
        labels{i} = sprintf('p%g', 100.0 * double(percentile(i)));
    end
end

function [work_avg, lift_avg] = compute_unthresholded_region_time_means(cfg, mask)
    ncfile = cfg.scenario.standard_ncfile;
    ensure_input_exists(ncfile);

    [~, ~, time] = read_standard_coordinates(ncfile);
    work3d = read_required_var(ncfile, 'work');
    lift3d = read_required_var(ncfile, 'lift');

    work_region = work3d(mask.lon_idx, mask.lat_idx, :);
    lift_region = lift3d(mask.lon_idx, mask.lat_idx, :);
    weight_3d = reshape(mask.lat_weights, 1, [], 1);

    work_num = sum(work_region .* weight_3d, [1, 2], 'omitnan');
    work_den = sum(double(~isnan(work_region)) .* weight_3d, [1, 2], 'omitnan');
    work_spatial_avg = squeeze(work_num ./ work_den);

    lift_num = sum(lift_region .* weight_3d, [1, 2], 'omitnan');
    lift_den = sum(double(~isnan(lift_region)) .* weight_3d, [1, 2], 'omitnan');
    lift_spatial_avg = squeeze(lift_num ./ lift_den);

    time_weights_days = cfg.scenario.time_weight_fn(time, ncfile);
    if iscell(time_weights_days)
        time_weights_days = time_weights_days{1};
    end

    work_avg = weighted_nanmean(work_spatial_avg, time_weights_days);
    lift_avg = weighted_nanmean(lift_spatial_avg, time_weights_days);
end

function [lat, lon, time] = read_standard_coordinates(ncfile)
    try
        lat = ncread(ncfile, 'lat');
        lon = ncread(ncfile, 'lon');
        time = ncread(ncfile, 'time');
    catch ME
        error('Error reading standard coordinates: %s', ME.message);
    end
end

function values = read_required_var(ncfile, varname)
    try
        values = ncread(ncfile, varname);
    catch ME
        error('Error reading %s variable: %s\nAvailable variables: %s', varname, ME.message, get_available_vars(ncfile));
    end
end

function stats = build_percentile_bin_stats(percentile, work_tail_avg, lift_tail_avg, work_total_avg, lift_total_avg)
    p = double(percentile(:));
    nperc = numel(p);
    if nperc < 1
        error('No percentile thresholds found in the thresholded file.');
    end
    if any(diff(p) <= 0)
        error('Percentile coordinate must be strictly increasing.');
    end
    if p(1) <= 0 || p(end) >= 1
        error('Percentile coordinate values must be strictly inside (0, 1).');
    end

    lower = [0.0; p];
    upper = [p; 1.0];
    mass = upper - lower;
    nbin = numel(mass);

    work_bin_avg = nan(nbin, 1);
    lift_bin_avg = nan(nbin, 1);

    tail_mass_first = 1.0 - p(1);
    work_bin_avg(1) = (work_total_avg - tail_mass_first * work_tail_avg(1)) / mass(1);
    lift_bin_avg(1) = (lift_total_avg - tail_mass_first * lift_tail_avg(1)) / mass(1);

    for ib = 2:nperc
        p_lo = p(ib - 1);
        p_hi = p(ib);
        tail_mass_lo = 1.0 - p_lo;
        tail_mass_hi = 1.0 - p_hi;
        work_bin_avg(ib) = (tail_mass_lo * work_tail_avg(ib - 1) - tail_mass_hi * work_tail_avg(ib)) / mass(ib);
        lift_bin_avg(ib) = (tail_mass_lo * lift_tail_avg(ib - 1) - tail_mass_hi * lift_tail_avg(ib)) / mass(ib);
    end

    work_bin_avg(end) = work_tail_avg(end);
    lift_bin_avg(end) = lift_tail_avg(end);

    ratio_bin = lift_bin_avg ./ work_bin_avg;
    ratio_bin(abs(work_bin_avg) < eps) = NaN;

    labels = cell(nbin, 1);
    labels{1} = sprintf('<P%s', format_percentile_value(100.0 * p(1)));
    for ib = 2:nperc
        labels{ib} = sprintf('P%s-P%s', format_percentile_value(100.0 * p(ib - 1)), format_percentile_value(100.0 * p(ib)));
    end
    labels{end} = sprintf('>P%s', format_percentile_value(100.0 * p(end)));

    stats = struct();
    stats.bin_lower = lower;
    stats.bin_upper = upper;
    stats.bin_mass = mass;
    stats.bin_center = 0.5 * (lower + upper);
    stats.work_bin_avg = work_bin_avg;
    stats.lift_bin_avg = lift_bin_avg;
    stats.ratio_bin = ratio_bin;
    stats.labels = labels;
end

function plot_percentile_bin_ratio(stats, region_title)
    nbin = numel(stats.ratio_bin);
    x = 1:nbin;
    y = stats.ratio_bin;

    figure('Position', [140, 140, 1050, 420]);
    set(gcf, 'Color', 'w');
    set(gcf, 'WindowStyle', 'docked');

    bar(x, y, 0.85, 'FaceColor', [0.20, 0.45, 0.70], 'EdgeColor', [0.05, 0.15, 0.25]);
    hold on;
    plot(x, y, 'k.-', 'LineWidth', 1.0, 'MarkerSize', 12);
    hold off;

    ax = gca;
    ax.XTick = x;
    ax.XTickLabel = stats.labels;
    ax.XTickLabelRotation = 35;
    ax.TickLabelInterpreter = 'none';
    xlim([0.5, nbin + 0.5]);

    xlabel('Percentile-rank bin');
    ylabel('Lift/Work ratio');
    title(sprintf('%s: Lift/Work Ratio by Percentile-Rank Bin', region_title), 'Interpreter', 'none');
    grid on;
    box on;
end

function print_percentile_bin_summary(stats)
    fprintf('\nPercentile-rank bin summary (surface-fraction aware):\n');
    for ib = 1:numel(stats.labels)
        fprintf('  %-15s  area=%.6f  work=%.6g  lift=%.6g  ratio=%.6g\n', ...
            stats.labels{ib}, stats.bin_mass(ib), stats.work_bin_avg(ib), stats.lift_bin_avg(ib), stats.ratio_bin(ib));
    end
end

function text_out = format_percentile_value(value_in)
    text_out = regexprep(sprintf('%.6f', value_in), '0+$', '');
    text_out = regexprep(text_out, '\.$', '');
end

function data_smoothed = running_mean_omitnan_cols(data_in, window_days)
    [ntime, ncol] = size(data_in);
    half_win = floor(window_days / 2);
    data_smoothed = nan(ntime, ncol);
    for icol = 1:ncol
        column = data_in(:, icol);
        for it = 1:ntime
            start_idx = max(1, it - half_win);
            stop_idx = min(ntime, it + half_win);
            segment = column(start_idx:stop_idx);
            valid = isfinite(segment);
            if any(valid)
                data_smoothed(it, icol) = mean(segment(valid));
            end
        end
    end
end