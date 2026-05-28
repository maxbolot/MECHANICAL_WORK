function out = run_control_work_lift_latband_bar_chart_panel_a_boxplot_analysis(cfg)
    cfg = apply_defaults(cfg);

    stats = compute_control_work_lift_latband_bar_chart_panel_a_boxplot_stats(cfg);
    stats = filter_stats_for_run_mode_local(stats, cfg.run_mode);
    if cfg.print_summary
        print_summary_local(stats);
    end
    if cfg.make_plots
        plot_panel_a_with_boxplots_local(stats, cfg);
    end

    out = stats;
end


function cfg = apply_defaults(cfg)
    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    if ~isfield(cfg, 'scenario') || isempty(cfg.scenario)
        cfg.scenario = scenario_control();
    end
    if ~isfield(cfg, 'threshold_percentile') || ~isfinite(cfg.threshold_percentile)
        cfg.threshold_percentile = 0.99999;
    end
    if ~isfield(cfg, 'make_plots')
        cfg.make_plots = true;
    end
    if ~isfield(cfg, 'print_summary')
        cfg.print_summary = true;
    end
    if ~isfield(cfg, 'window_style') || isempty(cfg.window_style)
        cfg.window_style = 'docked';
    end
    if ~isfield(cfg, 'run_mode') || isempty(cfg.run_mode)
        cfg.run_mode = 'tropics';
    end
    % ylim: empty means auto; supply a [lo hi] vector to override
    if ~isfield(cfg, 'ylim')
        cfg.ylim = [];
    end
end


function stats = filter_stats_for_run_mode_local(stats, run_mode)
    run_mode = string(run_mode);
    switch run_mode
        case "tropics"
            keep_name = 'tropics';
        case "midlatitudes"
            keep_name = 'midlatitudes';
        otherwise
            error('Unsupported run_mode "%s". Expected "tropics" or "midlatitudes".', run_mode);
    end

    keep_idx = find(strcmp({stats.regions.name}, keep_name));
    if isempty(keep_idx)
        error('Could not find region "%s" in stats.', keep_name);
    end

    stats.regions = stats.regions(keep_idx);
    stats.rows = stats.rows(keep_idx);
end


function print_summary_local(stats)
    fprintf('\n=== Control Work/Lift Lat-Band Panel A + Boxplot Summary ===\n');
    fprintf('Requested threshold percentile: %.8f\n', stats.threshold_percentile_requested);
    fprintf('Selected threshold percentile:  %.8f\n', stats.threshold_percentile_selected);
    fprintf('%-14s %12s %12s %14s %14s %14s %14s\n', ...
        'Region', 'thr work', 'thr lift', 'work out min', 'work out max', 'lift out min', 'lift out max');
    for ireg = 1:numel(stats.rows)
        row = stats.rows(ireg);
        [work_out_min, work_out_max] = outlier_minmax_local(row.threshold_work_samples);
        [lift_out_min, lift_out_max] = outlier_minmax_local(row.threshold_lift_samples);
        fprintf('%-14s %12.6g %12.6g %14s %14s %14s %14s\n', ...
            stats.regions(ireg).label, row.threshold_work, row.threshold_lift, ...
            format_value_local(work_out_min), format_value_local(work_out_max), ...
            format_value_local(lift_out_min), format_value_local(lift_out_max));
    end
end


function [out_min, out_max] = outlier_minmax_local(samples)
    valid = samples(isfinite(samples));
    if isempty(valid)
        out_min = NaN;
        out_max = NaN;
        return;
    end

    q1 = prctile(valid, 25);
    q3 = prctile(valid, 75);
    iqr_val = q3 - q1;
    lower_fence = q1 - 1.5 * iqr_val;
    upper_fence = q3 + 1.5 * iqr_val;
    outliers = valid(valid < lower_fence | valid > upper_fence);

    if isempty(outliers)
        out_min = NaN;
        out_max = NaN;
        return;
    end

    out_min = min(outliers);
    out_max = max(outliers);
end


function txt = format_value_local(v)
    if ~isfinite(v)
        txt = 'n/a';
        return;
    end
    txt = sprintf('%.6g', v);
end


function plot_panel_a_with_boxplots_local(stats, cfg)
    fig = figure('Units', 'inches', 'Position', [0, 0, 8.5, 4.8]);
    set(fig, 'Color', 'w');
    set(fig, 'WindowStyle', cfg.window_style);

    ax = axes(fig);

    row_labels = cell(numel(stats.rows), 1);
    panel_a_values = zeros(numel(stats.rows), 2);
    sample_cell = cell(numel(stats.rows), 2);
    for ireg = 1:numel(stats.rows)
        row = stats.rows(ireg);
        row_labels{ireg} = stats.regions(ireg).label;
        panel_a_values(ireg, :) = [row.threshold_work, row.threshold_lift];
        sample_cell{ireg, 1} = row.threshold_work_samples;
        sample_cell{ireg, 2} = row.threshold_lift_samples;
    end

    threshold_work_color = [0.82, 0.16, 0.16];
    threshold_lift_color = [0.10, 0.35, 0.80];

    [b1, x] = plot_grouped_bars_local(ax, panel_a_values, row_labels, 2, 0.16, 0.56, 2.10);
    b1(1).FaceColor = threshold_work_color;
    b1(2).FaceColor = threshold_lift_color;

    hold(ax, 'on');
    plot_boxplots_overlay_local(ax, sample_cell, x);
    hdist = plot(ax, nan, nan, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Distribution (boxplot)');
    hold(ax, 'off');

    % boxplot replaces XTick/XTickLabel — restore them
    xticks(ax, mean(x, 2));
    xticklabels(ax, row_labels);

    set(ax, 'YScale', 'linear');
    if ~isempty(cfg.ylim)
        ylim(ax, cfg.ylim);
    end
    ylabel(ax, 'W m^{-2}');
    title(ax, sprintf('a. P>%s work and lift by latitude family', percentile_label_local(stats.threshold_percentile_selected)), 'FontWeight', 'normal');
    set(ax, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', 'Box', 'off');

    legend(ax, [b1(1), b1(2), hdist], ...
        {sprintf('P>%s work', percentile_label_local(stats.threshold_percentile_selected)), ...
         sprintf('P>%s lift', percentile_label_local(stats.threshold_percentile_selected)), 'Distribution (boxplot)'}, ...
        'Location', 'northoutside', 'Orientation', 'horizontal');
    legend(ax, 'boxoff');

    sgtitle(sprintf('%s: control work/lift regional summary', stats.scenario.label), 'Interpreter', 'none');
end


function plot_boxplots_overlay_local(ax, sample_cell, x_positions)
    % Use consecutive group indices 1..N so boxplot's sorted-group-to-position
    % mapping is unambiguous: group k always maps to active_positions(k).
    [nclusters, nseries] = size(sample_cell);
    data = [];
    groups = [];
    active_positions = [];
    g = 0;

    for ic = 1:nclusters
        for is = 1:nseries
            samples = sample_cell{ic, is};
            samples = samples(isfinite(samples));
            if isempty(samples)
                continue;
            end
            g = g + 1;
            data = [data; samples(:)]; %#ok<AGROW>
            groups = [groups; g * ones(numel(samples), 1)]; %#ok<AGROW>
            active_positions = [active_positions; x_positions(ic, is)]; %#ok<AGROW>
        end
    end

    if isempty(data)
        return;
    end

    boxplot(ax, data, groups, ...
        'Positions', active_positions, ...
        'Widths', 0.10, ...
        'Colors', [0, 0, 0], ...
        'Symbol', 'k.');
end


function [bh, x] = plot_grouped_bars_local(ax, values, labels, nseries, bar_width, intra_cluster_spacing, cluster_spacing)
    axes(ax);
    nclusters = size(values, 1);
    if size(values, 2) ~= nseries
        error('Unexpected number of series: got %d, expected %d.', size(values, 2), nseries);
    end

    x = zeros(nclusters, nseries);
    half_span = ((nseries - 1) / 2) * intra_cluster_spacing;
    for ic = 1:nclusters
        center = 1 + (ic - 1) * cluster_spacing;
        offsets = ((1:nseries) - (nseries + 1) / 2) * intra_cluster_spacing;
        x(ic, :) = center + offsets;
    end

    hold(ax, 'on');
    bh = gobjects(1, nseries);
    for is = 1:nseries
        bh(is) = bar(ax, x(:, is), values(:, is), bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');
    end
    hold(ax, 'off');
    xlim(ax, [min(x(:)) - half_span - 0.35, max(x(:)) + half_span + 0.35]);
    xticks(ax, mean(x, 2));
    xticklabels(ax, labels);
end


function label = percentile_label_local(percentile_value)
    pct = 100.0 * percentile_value;
    pct_str = sprintf('%.10f', pct);
    pct_str = regexprep(pct_str, '0+$', '');
    pct_str = regexprep(pct_str, '\.$', '');
    label = ['P', strrep(pct_str, '.', '')];
end
