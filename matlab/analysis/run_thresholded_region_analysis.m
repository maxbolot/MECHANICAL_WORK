function out = run_thresholded_region_analysis(cfg)
    % Executes thresholded regional work/lift analysis and summary visualizations.
    ncfile = cfg.scenario.thresholded_ncfile;
    ensure_input_exists(ncfile);

    if ~isfield(cfg, 'smooth_before_plotting')
        cfg.smooth_before_plotting = false;
    end
    if ~isfield(cfg, 'smooth_window_days')
        cfg.smooth_window_days = 5;
    end
    if ~isfield(cfg, 'make_plots')
        cfg.make_plots = true;
    end
    if ~isfield(cfg, 'print_summary')
        cfg.print_summary = true;
    end

    fprintf('Reading NetCDF file: %s\n', ncfile);
    [lat, lon, time, percentile] = read_thresholded_coordinates(ncfile);
    fprintf('File dimensions: lat=%d, lon=%d, time=%d, percentile=%d\n', length(lat), length(lon), length(time), length(percentile));

    mask = build_region_mask(lat, lon, cfg.region);
    fprintf('%s\n', mask.selection_message);

    ts = compute_thresholded_region_spatial_timeseries(ncfile, mask);
    fprintf('Work variable reordered shape [lon lat time percentile]: %s\n', mat2str(ts.work4d_size));
    fprintf('Lift variable reordered shape [lon lat time percentile]: %s\n', mat2str(ts.lift4d_size));

    ntime = ts.ntime;
    nperc = ts.nperc;
    work_spatial_avg = ts.work_spatial_avg;
    lift_spatial_avg = ts.lift_spatial_avg;

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

    if cfg.print_summary
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

        print_percentile_bin_summary(bin_stats);
    end

    if cfg.make_plots
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
    end

    out = struct();
    out.percentile = percentile;
    out.percentile_labels = percentile_labels;
    out.time = time;
    out.missing_steps = missing_steps;
    out.work_spatial_avg = work_spatial_avg;
    out.lift_spatial_avg = lift_spatial_avg;
    out.work_avg = work_avg;
    out.lift_avg = lift_avg;
    out.ratio_avg = ratio_avg;
    out.bin_stats = bin_stats;
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
