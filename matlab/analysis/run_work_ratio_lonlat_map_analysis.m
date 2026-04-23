function out = run_work_ratio_lonlat_map_analysis(cfg)
    % Runs work/lift-ratio lon-lat analysis for selected scenario and time mode.
    cfg = apply_defaults(cfg);

    [ncfile, simulation_label, is_thresholded] = resolve_work_ratio_input_file(cfg.simulation);
    [lon, lat, time_axis] = read_coordinates(ncfile);
    time_dt = netcdf_time_to_datetime(time_axis, ncfile);

    [time_indices, time_label] = select_time_indices( ...
        time_dt, cfg.temporal_mode, cfg.instantaneous_date, cfg.range_start, cfg.range_end);

    if is_thresholded
        [work_data, selected_percentile] = read_lon_lat_time_from_var(ncfile, 'work', cfg.threshold_percentile);
        [lift_data, ~] = read_lon_lat_time_from_var(ncfile, 'lift', cfg.threshold_percentile);
        simulation_label = sprintf('%s (thresholded p%g)', simulation_label, 100.0 * selected_percentile);
    else
        work_data = read_required_var(ncfile, 'work');
        lift_data = read_required_var(ncfile, 'lift');
    end

    work_map = mean_over_time(work_data, time_indices);
    lift_map = mean_over_time(lift_data, time_indices);

    ratio_map = lift_map ./ work_map;
    ratio_map(abs(work_map) < eps) = NaN;

    if cfg.make_plots
        plot_maps(lon, lat, work_map, ratio_map, simulation_label, time_label);
    end

    out = struct();
    out.ncfile = ncfile;
    out.simulation = cfg.simulation;
    out.simulation_label = simulation_label;
    out.is_thresholded = is_thresholded;
    out.temporal_mode = cfg.temporal_mode;
    out.time_indices = time_indices;
    out.time_label = time_label;
    out.lon = lon;
    out.lat = lat;
    out.work_map = work_map;
    out.lift_map = lift_map;
    out.ratio_map = ratio_map;
end


function cfg = apply_defaults(cfg)
    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    if ~isfield(cfg, 'simulation')
        cfg.simulation = 'control';
    end
    if ~isfield(cfg, 'temporal_mode')
        cfg.temporal_mode = 'instantaneous';
    end
    if ~isfield(cfg, 'instantaneous_date')
        cfg.instantaneous_date = datetime(2021, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
    end
    if ~isfield(cfg, 'range_start')
        cfg.range_start = datetime(2020, 8, 1, 0, 0, 0, 'TimeZone', 'UTC');
    end
    if ~isfield(cfg, 'range_end')
        cfg.range_end = datetime(2020, 8, 30, 23, 59, 59, 'TimeZone', 'UTC');
    end
    if ~isfield(cfg, 'threshold_percentile')
        cfg.threshold_percentile = 0.5;
    end
    if ~isfield(cfg, 'make_plots')
        cfg.make_plots = true;
    end
end


function [lon, lat, time_axis] = read_coordinates(ncfile)
    lon = ncread(ncfile, 'lon');
    lat = ncread(ncfile, 'lat');
    time_axis = ncread(ncfile, 'time');
end


function field_map = mean_over_time(field_data, time_indices)
    if ndims(field_data) ~= 3
        error('Expected field_data dimensions [lon, lat, time], got %s', mat2str(size(field_data)));
    end
    field_map = mean(field_data(:, :, time_indices), 3, 'omitnan');
end


function plot_maps(lon, lat, work_map, ratio_map, simulation_label, time_label)
    figure('Position', [100, 100, 1300, 500]);
    set(gcf, 'Color', 'w');
    set(gcf, 'WindowStyle', 'docked');

    t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    use_mapping = has_mapping_toolbox();
    [x_grid, y_grid] = ndgrid(double(lon(:)), double(lat(:)));

    ax1 = nexttile;
    plot_panel_map(ax1, x_grid, y_grid, work_map, use_mapping);
    title(ax1, sprintf('a. Work (%s)', time_label), 'Interpreter', 'none');
    cb1 = colorbar(ax1, 'eastoutside');
    cb1.Label.String = 'W m^{-2}';
    colormap(ax1, parula);
    clim(ax1, [0, 50]);

    ax2 = nexttile;
    plot_panel_map(ax2, x_grid, y_grid, ratio_map, use_mapping);
    title(ax2, sprintf('b. Lift/Work (%s)', time_label), 'Interpreter', 'none');
    cb2 = colorbar(ax2, 'eastoutside');
    cb2.Label.String = 'unitless';
    colormap(ax2, parula);
    clim(ax2, [0, 1]);

    sgtitle(t, sprintf('%s: Work and Lift/Work in Lon-Lat Coordinates', simulation_label), 'Interpreter', 'none');
end


function plot_panel_map(ax, x_grid, y_grid, field_map, use_mapping)
    if use_mapping
        axes(ax);
        axesm('mercator', 'Frame', 'on', 'Grid', 'on', 'MapLatLimit', [-60 60], 'MapLonLimit', [-180 180], ...
            'meridianlabel', 'on', 'parallellabel', 'on', 'labelformat', 'signed', ...
            'mlabellocation', 60, 'plinelocation', 20, 'mlinelocation', 45, ...
            'MLabelParallel', 'south', 'FontSize', 10);
        geoshow(y_grid, x_grid, field_map, 'DisplayType', 'texturemap');
        shading flat;
        try
            geoshow('landareas.shp', 'FaceColor', 'none');
        catch
        end
        axis off;
        set(gca, 'Color', 'none');
        set(gca, 'Layer', 'top');
        return;
    end

    imagesc(ax, x_grid(:, 1), y_grid(1, :), field_map.');
    set(ax, 'YDir', 'normal');
    xlim(ax, [-180 180]);
    ylim(ax, [-60 60]);
    xlabel(ax, 'Longitude');
    ylabel(ax, 'Latitude');
end


function tf = has_mapping_toolbox()
    tf = exist('axesm', 'file') == 2 && exist('geoshow', 'file') == 2;
end
