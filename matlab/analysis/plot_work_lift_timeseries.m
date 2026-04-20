function plot_work_lift_timeseries(plot_cfg)
    figure('Position', [100, 100, 1100, 700]);
    set(gcf, 'color', 'w');
    set(gcf, 'WindowStyle', 'docked');

    if isfield(plot_cfg, 'legend_labels') && ~isempty(plot_cfg.legend_labels)
        line_colors = lines(numel(plot_cfg.legend_labels));
    else
        line_colors = [];
    end

    subplot(3, 1, 1);
    plot_series(plot_cfg.plot_time, plot_cfg.work_plot, plot_cfg.plot_time_label, 'Work', ...
        sprintf('Mechanical Work - %s Average Time Series%s', plot_cfg.title_prefix, plot_cfg.title_suffix), line_colors);

    subplot(3, 1, 2);
    plot_series(plot_cfg.plot_time, plot_cfg.lift_plot, plot_cfg.plot_time_label, 'Lift', ...
        sprintf('Lift Work - %s Average Time Series%s', plot_cfg.title_prefix, plot_cfg.title_suffix), line_colors);

    subplot(3, 1, 3);
    plot_series(plot_cfg.plot_time, plot_cfg.ratio_plot, plot_cfg.plot_time_label, 'Lift/Work', ...
        sprintf('Lift/Work Ratio - %s Average Time Series%s', plot_cfg.title_prefix, plot_cfg.title_suffix), line_colors);
    ylim([0 1]);

    if ~isempty(line_colors)
        lgd = legend(plot_cfg.legend_labels, 'Location', 'eastoutside');
        lgd.Title.String = 'Threshold';
    end
end

function plot_series(plot_time, values, plot_time_label, y_label, plot_title, line_colors)
    if isvector(values)
        plot(plot_time, values, 'LineWidth', 1.5);
    else
        hold on;
        for i = 1:size(values, 2)
            plot(plot_time, values(:, i), 'LineWidth', 1.2, 'Color', line_colors(i, :));
        end
        hold off;
    end
    xlabel(plot_time_label);
    ylabel(y_label);
    title(plot_title);
    grid on;
end