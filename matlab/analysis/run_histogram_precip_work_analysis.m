function out = run_histogram_precip_work_analysis(cfg)
    % Orchestrates tropical precipitation-histogram diagnostics and optional plotting.
    cfg = apply_defaults(cfg);

    ensure_input_exists(cfg.ncfile);
    fprintf('Reading NetCDF file: %s\n', cfg.ncfile);

    hist_data = read_tropical_histogram_diagnostics(cfg.ncfile);

    [time_weights_days, missing_steps] = cfg.time_weight_fn(hist_data.time, cfg.ncfile);
    metrics = compute_histogram_precip_work_metrics(hist_data, time_weights_days);

    if cfg.print_summary
        fprintf('\n=== RESULTS ===\n');
        fprintf('Simulation: %s\n', cfg.simulation_name);
        fprintf('Time weighting: uniform 1-day cadence\n');
        fprintf('Detected missing timesteps: %d\n', missing_steps);
        fprintf('Number of time steps: %d\n', hist_data.ntime);
    end

    if cfg.make_plots
        plot_histogram_outputs(hist_data, metrics, cfg);
    end

    out = struct();
    out.hist_data = hist_data;
    out.metrics = metrics;
    out.missing_steps = missing_steps;
    out.simulation_name = cfg.simulation_name;
end


function cfg = apply_defaults(cfg)
    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    if ~isfield(cfg, 'simulation_name')
        cfg.simulation_name = 'control';
    end
    if ~isfield(cfg, 'ncfile')
        cfg.ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms/hist_2020010300_2022011200.nc';
    end
    if ~isfield(cfg, 'time_weight_fn')
        cfg.time_weight_fn = @compute_time_weights_control;
    end
    if ~isfield(cfg, 'xmax')
        cfg.xmax = 1.5e-2;
    end
    if ~isfield(cfg, 'title_prefix')
        cfg.title_prefix = 'Control';
    end
    if ~isfield(cfg, 'make_plots')
        cfg.make_plots = true;
    end
    if ~isfield(cfg, 'print_summary')
        cfg.print_summary = true;
    end
end


function plot_histogram_outputs(hist_data, metrics, cfg)
    fig = figure('Position', [120, 120, 1000, 700]);
    set(fig, 'color', 'w');
    set(fig, 'WindowStyle', 'docked');

    subplot(2,1,1);
    pr_freq_plot = metrics.pr_freq_timeavg;
    pr_freq_plot(pr_freq_plot <= 0) = NaN;
    plot(hist_data.pr_centers, pr_freq_plot, 'k-', 'LineWidth', 1.8);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([min(hist_data.pr_centers), cfg.xmax]);
    xlabel('Precipitation bin center');
    ylabel('Frequency');
    title(sprintf('%s: Time-Mean Area-Weighted Precipitation Frequency', cfg.title_prefix));
    grid on;

    subplot(2,1,2);
    lift_plot = metrics.lift_by_pr_timeavg;
    work_plot = metrics.work_by_pr_timeavg;
    lift_plot(lift_plot <= 0) = NaN;
    work_plot(work_plot <= 0) = NaN;
    ratio_plot = metrics.ratio_by_pr_timeavg;
    ratio_plot(~isfinite(ratio_plot) | ratio_plot <= 0) = NaN;

    yyaxis left;
    h1 = plot(hist_data.pr_centers, lift_plot, 'r-', 'LineWidth', 1.8); hold on;
    h2 = plot(hist_data.pr_centers, work_plot, 'b-', 'LineWidth', 1.8);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    ylabel('Power density (W/m^2)');

    yyaxis right;
    h3 = plot(hist_data.pr_centers, ratio_plot, 'k--', 'LineWidth', 1.5);
    ylabel('Lift/Work ratio');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

    xlim([min(hist_data.pr_centers), cfg.xmax]);
    xlabel('Precipitation bin center');
    title(sprintf('%s: Time-Mean Lift and Work in Precipitation Bins', cfg.title_prefix));
    legend([h1 h2 h3], {'Lift', 'Work', 'Lift/Work'}, 'Location', 'best');
    grid on;
    hold off;
end
