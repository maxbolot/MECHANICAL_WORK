function out = run_work_lift_all_regions_analysis(cfg)
    cfg = apply_defaults(cfg);

    ensure_input_exists(cfg.control_file);
    ensure_input_exists(cfg.warming_file);

    fprintf('Reading control simulation: %s\n', cfg.control_file);
    control = process_work_lift_regions_file(cfg.control_file, cfg.regions, 'control');

    fprintf('Reading warming simulation: %s\n', cfg.warming_file);
    warming = process_work_lift_regions_file(cfg.warming_file, cfg.regions, 'warming');

    % Build stacked-bar input [lift, ke] in legacy row order:
    % [control zone1; warming zone1; control zone2; warming zone2; ...]
    Y = zeros(8, 2);
    for i = 1:4
        Y(2*i - 1, :) = [control.lift(i), control.ke(i)];
        Y(2*i, :) = [warming.lift(i), warming.ke(i)];
    end

    % Add spacer row after each region pair to mimic legacy axis 1 layout.
    newY = zeros(12, 2);
    newY([1 2 4 5 7 8 10 11], :) = Y;

    if cfg.make_plots
        plot_legacy_bars(newY);
    end

    if cfg.print_summary
        fprintf('\n=== Spatial + Time Averaged Results ===\n');
        fprintf('%-22s %-12s %12s %12s %12s %12s %10s\n', 'Region', 'Simulation', 'work', 'lift', 'ke', 'lift/work', 'missing');
        for i = 1:4
            fprintf('%-22s %-12s %12.6f %12.6f %12.6f %12.6f %10d\n', cfg.regions{i}.name, 'control', ...
                control.work(i), control.lift(i), control.ke(i), control.ratio(i), control.missing_steps(i));
            fprintf('%-22s %-12s %12.6f %12.6f %12.6f %12.6f %10d\n', cfg.regions{i}.name, 'warming', ...
                warming.work(i), warming.lift(i), warming.ke(i), warming.ratio(i), warming.missing_steps(i));
        end
    end

    out = struct();
    out.control = control;
    out.warming = warming;
    out.newY = newY;
    out.regions = cfg.regions;
end


function cfg = apply_defaults(cfg)
    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    if ~isfield(cfg, 'control_file')
        cfg.control_file = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180/work_2020010300_2022011200.nc';
    end
    if ~isfield(cfg, 'warming_file')
        cfg.warming_file = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/work_2020010300_2022011800.nc';
    end
    if ~isfield(cfg, 'regions')
        cfg.regions = { ...
            struct('name', 'global', 'latmin', -90, 'latmax', 90, 'lonmin', NaN, 'lonmax', NaN), ...
            struct('name', 'tropics', 'latmin', -30, 'latmax', 30, 'lonmin', NaN, 'lonmax', NaN), ...
            struct('name', 'maritime continent', 'latmin', -15, 'latmax', 15, 'lonmin', 90, 'lonmax', 150), ...
            struct('name', 'north. midlatitudes', 'latmin', 30, 'latmax', 60, 'lonmin', NaN, 'lonmax', NaN) ...
        };
    end
    if ~isfield(cfg, 'make_plots')
        cfg.make_plots = true;
    end
    if ~isfield(cfg, 'print_summary')
        cfg.print_summary = true;
    end
end


function plot_legacy_bars(newY)
    fig = figure('units', 'inch', 'position', [0, 0, 8, 3.8]);
    set(fig, 'color', 'w');
    set(fig, 'WindowStyle', 'docked');

    ax1 = axes(fig);
    set(ax1, 'FontSize', 13);
    b1 = bar(newY, 'stacked', 'FaceColor', 'flat');
    set(gca, 'XLim', [0 12], 'XTick', [1.5 4.5 7.5 10.5], ...
        'XTickLabel', {'global', 'tropics', 'maritime continent', 'north. midlatitudes'});
    ylabel('W m^{-2}');
    title('Legacy Axis 1: work = lift + ke');

    control_dark = [0.10, 0.35, 0.80];
    control_light = [0.65, 0.80, 0.98];
    warming_dark = [0.85, 0.20, 0.20];
    warming_light = [0.98, 0.72, 0.72];
    white_bar = [1.0, 1.0, 1.0];

    c1 = [control_dark; warming_dark; white_bar; control_dark; warming_dark; white_bar; ...
          control_dark; warming_dark; white_bar; control_dark; warming_dark; white_bar];
    c2 = [control_light; warming_light; white_bar; control_light; warming_light; white_bar; ...
          control_light; warming_light; white_bar; control_light; warming_light; white_bar];

    b1(1).CData = c1;
    b1(2).CData = c2;

    legend({'lift', 'ke = work - lift'}, 'Location', 'northwest');
    legend('boxoff');
    grid on;
end
