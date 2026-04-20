function plot_percentile_bin_ratio_comparison(stats_ctrl, stats_warm, region_title, label_ctrl, label_warm)
% PLOT_PERCENTILE_BIN_RATIO_COMPARISON  Superimposed grouped-bar chart of
%   Lift/Work ratio for two scenarios (control in blue, warming in red).
%
%   plot_percentile_bin_ratio_comparison(stats_ctrl, stats_warm, region_title)
%   plot_percentile_bin_ratio_comparison(stats_ctrl, stats_warm, region_title, label_ctrl, label_warm)
%
%   stats_ctrl / stats_warm : structs returned by compute_thresholded_bin_stats
%   region_title            : string used in the figure title
%   label_ctrl              : (optional) legend label for control, default 'Control'
%   label_warm              : (optional) legend label for warming, default '+4K'

    if nargin < 4 || isempty(label_ctrl)
        label_ctrl = 'Control';
    end
    if nargin < 5 || isempty(label_warm)
        label_warm = '+4K';
    end

    nbin   = numel(stats_ctrl.ratio_bin);
    x      = (1:nbin)';
    y_ctrl = stats_ctrl.ratio_bin;
    y_warm = stats_warm.ratio_bin;

    color_ctrl = [0.20, 0.45, 0.70];   % blue
    color_warm = [0.80, 0.20, 0.20];   % red
    edge_ctrl  = [0.05, 0.15, 0.25];
    edge_warm  = [0.35, 0.05, 0.05];

    bar_width = 0.75;
    alpha_ctrl = 0.75;
    alpha_warm = 1.00;

    figure('Position', [140, 140, 1050, 460]);
    set(gcf, 'Color', 'w');
    set(gcf, 'WindowStyle', 'docked');

    hold on;

    % Draw overlapping bars — warming first (behind), then control on top
    for ib = 1:nbin
        if isfinite(y_warm(ib))
            b = bar(x(ib), y_warm(ib), bar_width, ...
                'FaceColor', color_warm, 'EdgeColor', edge_warm);
            b.FaceAlpha = alpha_warm;
        end
        if isfinite(y_ctrl(ib))
            b = bar(x(ib), y_ctrl(ib), bar_width, ...
                'FaceColor', color_ctrl, 'EdgeColor', edge_ctrl);
            b.FaceAlpha = alpha_ctrl;
        end
    end

    % Overlay dot-line for each scenario
    h_ctrl = plot(x, y_ctrl, '.-', 'Color', color_ctrl * 0.7, ...
        'LineWidth', 1.2, 'MarkerSize', 12);
    h_warm = plot(x, y_warm, '.-', 'Color', color_warm * 0.7, ...
        'LineWidth', 1.2, 'MarkerSize', 12);

    hold off;

    ax = gca;
    ax.XTick       = x;
    ax.XTickLabel  = stats_ctrl.labels;
    ax.XTickLabelRotation = 35;
    ax.TickLabelInterpreter = 'none';
    xlim([0.5, nbin + 0.5]);

    legend([h_ctrl, h_warm], {label_ctrl, label_warm}, 'Location', 'northwest');
    xlabel('Percentile-rank bin');
    ylabel('Lift/Work ratio');
    title(sprintf('%s: Lift/Work Ratio by Percentile-Rank Bin (%s vs %s)', ...
        region_title, label_ctrl, label_warm), 'Interpreter', 'none');
    grid on;
    box on;
end
