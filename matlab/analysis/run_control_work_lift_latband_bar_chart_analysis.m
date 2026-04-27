function out = run_control_work_lift_latband_bar_chart_analysis(cfg)
    cfg = apply_defaults(cfg);

    stats = compute_control_work_lift_latband_bar_chart_stats(cfg);
    if cfg.print_summary
        print_summary_local(stats);
    end
    if cfg.make_plots
        plot_bar_chart_local(stats, cfg);
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
end


function print_summary_local(stats)
    fprintf('\n=== Control Work/Lift Lat-Band Bar Chart Summary ===\n');
    fprintf('Requested threshold percentile: %.8f\n', stats.threshold_percentile_requested);
    fprintf('Selected threshold percentile:  %.8f\n', stats.threshold_percentile_selected);
    fprintf('%-14s %12s %12s %12s %12s %12s %12s\n', ...
        'Region', 'mean work', 'mean lift', 'thr work', 'thr lift', 'mean L/W', 'thr L/W');
    for ireg = 1:numel(stats.rows)
        row = stats.rows(ireg);
        fprintf('%-14s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n', ...
            stats.regions(ireg).label, row.mean_work, row.mean_lift, ...
            row.threshold_work, row.threshold_lift, row.mean_ratio, row.threshold_ratio);
    end
end


function plot_bar_chart_local(stats, cfg)
    fig = figure('Units', 'inches', 'Position', [0, 0, 10.5, 4.8]);
    set(fig, 'Color', 'w');
    set(fig, 'WindowStyle', cfg.window_style);

    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    row_labels = cell(numel(stats.rows), 1);
    panel_a_values = zeros(numel(stats.rows), 4);
    panel_b_values = zeros(numel(stats.rows), 2);
    for ireg = 1:numel(stats.rows)
        row = stats.rows(ireg);
        row_labels{ireg} = stats.regions(ireg).label;
        panel_a_values(ireg, :) = [row.mean_work, row.mean_lift, row.threshold_work, row.threshold_lift];
        panel_b_values(ireg, :) = [row.mean_ratio, row.threshold_ratio];
    end

    mean_work_color = [0.98, 0.78, 0.78];
    mean_lift_color = [0.68, 0.82, 0.98];
    threshold_work_color = [0.82, 0.16, 0.16];
    threshold_lift_color = [0.10, 0.35, 0.80];

    ax1 = nexttile(1);
    b1 = plot_grouped_bars_local(ax1, panel_a_values, row_labels, 4, 0.17, 0.40, 2.10);
    set(ax1, 'YScale', 'log');
    ylim(ax1, [0.5, 1.0e2]);
    ylabel(ax1, 'W m^{-2}');
    title(ax1, 'a. Work and lift by latitude family', 'FontWeight', 'normal');
    b1(1).FaceColor = mean_work_color;
    b1(2).FaceColor = mean_lift_color;
    b1(3).FaceColor = threshold_work_color;
    b1(4).FaceColor = threshold_lift_color;
    legend(ax1, {'Mean work', 'Mean lift', sprintf('P>%s work', percentile_label_local(stats.threshold_percentile_selected)), sprintf('P>%s lift', percentile_label_local(stats.threshold_percentile_selected))}, ...
        'Location', 'northoutside', 'Orientation', 'horizontal');
    legend(ax1, 'boxoff');
    set(ax1, 'XGrid', 'off', 'YGrid', 'off');

    ax2 = nexttile(2);
    b2 = plot_grouped_bars_local(ax2, panel_b_values, row_labels, 2, 0.22, 0.48, 1.55);
    ylabel(ax2, 'Lift / Work');
    ylim(ax2, [0, 1.05]);
    title(ax2, 'b. Lift/work ratio by latitude family', 'FontWeight', 'normal');
    b2(1).FaceColor = [0.82, 0.82, 0.82];
    b2(2).FaceColor = [0.35, 0.35, 0.35];
    yline(ax2, 1.0, 'k--', 'LineWidth', 1.5, 'Label', 'Theoretical Maximum (W_{lift} = W_{tot})');
    legend(ax2, {'Mean', sprintf('P>%s', percentile_label_local(stats.threshold_percentile_selected))}, 'Location', 'northwest');
    legend(ax2, 'boxoff');
    set(ax2, 'XGrid', 'off', 'YGrid', 'off');

    sgtitle(t, sprintf('%s: control work/lift regional summary', stats.scenario.label), 'Interpreter', 'none');
end


function bh = plot_grouped_bars_local(ax, values, labels, nseries, bar_width, intra_cluster_spacing, cluster_spacing)
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