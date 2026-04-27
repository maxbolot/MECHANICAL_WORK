function out = run_tropics_precip_histogram_time_mean_analysis(cfg)
    % Computes and plots simulation-mean tropical precipitation histogram from hist_area.
    cfg = apply_defaults(cfg);

    if ~isfile(cfg.ncfile)
        error('Histogram file not found: %s', cfg.ncfile);
    end

    pr_edges = double(ncread(cfg.ncfile, 'pr_edges'));
    hist_area_raw = double(ncread(cfg.ncfile, 'hist_area'));
    time_axis = double(ncread(cfg.ncfile, 'time'));

    nbin_pr = numel(pr_edges) - 1;
    if nbin_pr <= 0
        error('Invalid pr_edges in %s (need at least 2 edges).', cfg.ncfile);
    end

    hist_area_time_bin = orient_hist_area_time_bin(hist_area_raw, nbin_pr);

    [time_weights_days, missing_steps] = cfg.time_weight_fn(time_axis, cfg.ncfile);
    mean_hist_area = weighted_time_mean_rows_local(hist_area_time_bin, time_weights_days);

    total_area = sum(mean_hist_area, 'omitnan');
    area_fraction = mean_hist_area ./ total_area;
    if ~(isfinite(total_area) && (total_area > 0))
        area_fraction(:) = NaN;
    end

    pr_bin_centers = 0.5 * (pr_edges(1:end-1) + pr_edges(2:end));

    if cfg.make_plot
        plot_precip_histograms_local(pr_bin_centers, mean_hist_area, area_fraction, cfg);
    end

    if cfg.print_summary
        fprintf('Loaded %s\n', cfg.ncfile);
        fprintf('Time samples: %d, precipitation bins: %d, missing daily steps: %d\n', ...
            size(hist_area_time_bin, 1), size(hist_area_time_bin, 2), missing_steps);
    end

    out = struct();
    out.pr_edges = pr_edges;
    out.pr_bin_centers = pr_bin_centers;
    out.mean_hist_area = mean_hist_area;
    out.area_fraction = area_fraction;
    out.missing_steps = missing_steps;
end


function cfg = apply_defaults(cfg)
    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    if ~isfield(cfg, 'ncfile')
        cfg.ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms/hist_2020010300_2022011200.nc';
    end
    if ~isfield(cfg, 'time_weight_fn')
        cfg.time_weight_fn = @compute_time_weights_control;
    end
    if ~isfield(cfg, 'title_prefix')
        cfg.title_prefix = 'Control';
    end
    if ~isfield(cfg, 'xmax')
        cfg.xmax = [];
    end
    if ~isfield(cfg, 'xlim')
        cfg.xlim = [];
    end
    if ~isfield(cfg, 'ylim_hist_area')
        cfg.ylim_hist_area = [];
    end
    if ~isfield(cfg, 'ylim_area_fraction')
        cfg.ylim_area_fraction = [];
    end
    if ~isfield(cfg, 'make_plot')
        cfg.make_plot = true;
    end
    if ~isfield(cfg, 'print_summary')
        cfg.print_summary = true;
    end
end


function hist_area_time_bin = orient_hist_area_time_bin(hist_area_raw, nbin_pr)
    sz = size(hist_area_raw);
    if numel(sz) ~= 2
        error('hist_area must be 2D [time, nbin_pr] (or transposed). Got size %s.', mat2str(sz));
    end

    if sz(2) == nbin_pr
        hist_area_time_bin = hist_area_raw;
    elseif sz(1) == nbin_pr
        hist_area_time_bin = hist_area_raw.';
    else
        error(['Could not locate precipitation-bin dimension in hist_area. ', ...
            'Expected one dimension to be nbin_pr=%d, got size %s.'], nbin_pr, mat2str(sz));
    end
end


function mean_vals = weighted_time_mean_rows_local(values_time_bin, time_weights_days)
    if size(values_time_bin, 1) ~= numel(time_weights_days)
        error(['Time dimension mismatch for hist_area and weights: ', ...
            'hist_area has %d rows, weights has %d elements.'], ...
            size(values_time_bin, 1), numel(time_weights_days));
    end

    weights = double(time_weights_days(:));
    valid = isfinite(values_time_bin);

    weighted_num = sum(values_time_bin .* weights, 1, 'omitnan');
    weighted_den = sum(double(valid) .* weights, 1, 'omitnan');

    mean_vals = weighted_num ./ weighted_den;
    mean_vals(weighted_den <= 0) = NaN;
end


function plot_precip_histograms_local(pr_bin_centers, mean_hist_area, area_fraction, cfg)
    fig = figure('Color', 'w', 'Position', [100, 100, 950, 700]);
    set(fig, 'WindowStyle', 'docked');
    t = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(t, 1);
    plot(ax1, pr_bin_centers, mean_hist_area, 'LineWidth', 1.5, 'Color', [0.0, 0.35, 0.75]);
    grid(ax1, 'on');
    xlabel(ax1, 'Precipitation bin center');
    ylabel(ax1, 'Mean hist_area');
    title(ax1, 'Tropics: time-mean precipitation histogram (area per bin)');

    if all(pr_bin_centers > 0)
        set(ax1, 'XScale', 'log');
    end
    x_limits = resolve_x_limits_local(pr_bin_centers, cfg);
    if ~isempty(x_limits)
        xlim(ax1, x_limits);
    end
    if ~isempty(cfg.ylim_hist_area)
        ylim(ax1, validate_limits_local(cfg.ylim_hist_area, 'ylim_hist_area'));
    end

    ax2 = nexttile(t, 2);
    plot(ax2, pr_bin_centers, area_fraction, 'LineWidth', 1.5, 'Color', [0.8, 0.25, 0.1]);
    grid(ax2, 'on');
    xlabel(ax2, 'Precipitation bin center');
    ylabel(ax2, 'Area fraction');
    title(ax2, 'Tropics: time-mean normalized precipitation histogram');

    if all(pr_bin_centers > 0)
        set(ax2, 'XScale', 'log');
    end
    if ~isempty(x_limits)
        xlim(ax2, x_limits);
    end
    if ~isempty(cfg.ylim_area_fraction)
        ylim(ax2, validate_limits_local(cfg.ylim_area_fraction, 'ylim_area_fraction'));
    end

    sgtitle(t, {sprintf('%s simulation', cfg.title_prefix), 'Tropics precipitation histogram (simulation mean)'}, ...
        'Interpreter', 'none');
end


function x_limits = resolve_x_limits_local(pr_bin_centers, cfg)
    if ~isempty(cfg.xlim)
        x_limits = validate_limits_local(cfg.xlim, 'xlim');
        return;
    end

    if ~isempty(cfg.xmax) && isfinite(cfg.xmax)
        x_limits = [min(pr_bin_centers), cfg.xmax];
        return;
    end

    x_limits = [];
end


function limits = validate_limits_local(values, field_name)
    if ~(isnumeric(values) && numel(values) == 2 && all(isfinite(values(:))))
        error('%s must be a finite 2-element numeric vector [min max].', field_name);
    end

    limits = double(values(:)).';
    if limits(1) >= limits(2)
        error('%s must satisfy min < max.', field_name);
    end
end
