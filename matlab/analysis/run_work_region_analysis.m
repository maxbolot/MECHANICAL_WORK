function run_work_region_analysis(cfg)
    % Computes unthresholded regional work/lift means and time-series plots.
    ncfile = cfg.scenario.standard_ncfile;
    ensure_input_exists(ncfile);

    fprintf('Reading NetCDF file: %s\n', ncfile);
    [lat, lon, time] = read_standard_coordinates(ncfile);
    fprintf('File dimensions: lat=%d, lon=%d, time=%d\n', length(lat), length(lon), length(time));

    mask = build_region_mask(lat, lon, cfg.region);
    fprintf('%s\n', mask.selection_message);

    work = read_required_var(ncfile, 'work');
    fprintf('Work variable shape: %s\n', mat2str(size(work)));
    lift = read_required_var(ncfile, 'lift');
    fprintf('Lift variable shape: %s\n', mat2str(size(lift)));

    work_region = work(mask.lon_idx, mask.lat_idx, :);
    lift_region = lift(mask.lon_idx, mask.lat_idx, :);
    weight_3d = reshape(mask.lat_weights, 1, [], 1);

    work_num = sum(work_region .* weight_3d, [1, 2], 'omitnan');
    work_den = sum(double(~isnan(work_region)) .* weight_3d, [1, 2], 'omitnan');
    work_spatial_avg = squeeze(work_num ./ work_den);

    lift_num = sum(lift_region .* weight_3d, [1, 2], 'omitnan');
    lift_den = sum(double(~isnan(lift_region)) .* weight_3d, [1, 2], 'omitnan');
    lift_spatial_avg = squeeze(lift_num ./ lift_den);

    [time_weights_days, missing_steps] = cfg.scenario.time_weight_fn(time, ncfile);
    work_avg = weighted_nanmean(work_spatial_avg, time_weights_days);
    lift_avg = weighted_nanmean(lift_spatial_avg, time_weights_days);

    ratio_spatial_avg = lift_spatial_avg ./ work_spatial_avg;
    ratio_spatial_avg(abs(work_spatial_avg) < eps) = NaN;
    ratio_avg = lift_avg / work_avg;
    if abs(work_avg) < eps
        ratio_avg = NaN;
    end

    [plot_time, plot_time_label] = build_time_axis(ncfile, time);

    summarize_work_lift(struct( ...
        'region_label', cfg.region.result_label, ...
        'average_label', cfg.region.average_label, ...
        'missing_steps', missing_steps, ...
        'ntime', numel(time), ...
        'work_avg', work_avg, ...
        'lift_avg', lift_avg, ...
        'ratio_avg', ratio_avg));

    empty_legend = cell(0, 1);
    plot_work_lift_timeseries(struct( ...
        'plot_time', plot_time, ...
        'plot_time_label', plot_time_label, ...
        'work_plot', work_spatial_avg, ...
        'lift_plot', lift_spatial_avg, ...
        'ratio_plot', ratio_spatial_avg, ...
        'title_prefix', cfg.region.title_label, ...
        'title_suffix', '', ...
        'legend_labels', {empty_legend}));
end
