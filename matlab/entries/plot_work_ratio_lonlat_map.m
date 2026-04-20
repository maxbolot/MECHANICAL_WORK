if ~exist('matlab_base_dir', 'var')
    this_file = mfilename('fullpath');
    script_dir = fileparts(this_file);
    matlab_base_dir = fileparts(script_dir);
end
addpath(fullfile(matlab_base_dir, 'lib'));
addpath(fullfile(matlab_base_dir, 'presets'));

% User switches
simulation = 'control';
% Supported: control | warming | control_prate_thresholded | warming_prate_thresholded

temporal_mode = 'instantaneous';
% Supported: instantaneous | timerange

% Preset times for display
instantaneous_date = datetime(2021, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
range_start = datetime(2020, 8, 1, 0, 0, 0, 'TimeZone', 'UTC');
range_end = datetime(2020, 8, 30, 23, 59, 59, 'TimeZone', 'UTC');

% Used only for thresholded simulations.
threshold_percentile = 0.5;

[ncfile, simulation_label, is_thresholded] = resolve_input_file(simulation);
[lon, lat, time_axis] = read_coordinates(ncfile);
time_dt = netcdf_time_to_datetime(ncfile, time_axis);

[time_indices, time_label] = select_time_indices(time_dt, temporal_mode, instantaneous_date, range_start, range_end);

if is_thresholded
    [work_data, selected_percentile] = read_thresholded_field_at_percentile(ncfile, 'work', threshold_percentile);
    [lift_data, ~] = read_thresholded_field_at_percentile(ncfile, 'lift', threshold_percentile);
    simulation_label = sprintf('%s (thresholded p%g)', simulation_label, 100.0 * selected_percentile);
else
    work_data = ncread(ncfile, 'work');
    lift_data = ncread(ncfile, 'lift');
end

work_map = mean_over_time(work_data, time_indices);
lift_map = mean_over_time(lift_data, time_indices);

ratio_map = lift_map ./ work_map;
ratio_map(abs(work_map) < eps) = NaN;

plot_maps(lon, lat, work_map, ratio_map, simulation_label, time_label);

function [ncfile, simulation_label, is_thresholded] = resolve_input_file(simulation)
    switch lower(string(simulation))
        case "control"
            sc = scenario_control();
            ncfile = sc.standard_ncfile;
            simulation_label = 'Control';
            is_thresholded = false;
        case "warming"
            sc = scenario_plus4k();
            ncfile = sc.standard_ncfile;
            simulation_label = 'Warming (+4K)';
            is_thresholded = false;
        case "control_prate_thresholded"
            sc = scenario_control();
            ncfile = sc.thresholded_ncfile;
            simulation_label = 'Control';
            is_thresholded = true;
        case "warming_prate_thresholded"
            sc = scenario_plus4k();
            ncfile = sc.thresholded_ncfile;
            simulation_label = 'Warming (+4K)';
            is_thresholded = true;
        otherwise
            error('Unsupported simulation: %s', simulation);
    end

    if ~isfile(ncfile)
        error('File not found: %s', ncfile);
    end
end

function [lon, lat, time_axis] = read_coordinates(ncfile)
    lon = ncread(ncfile, 'lon');
    lat = ncread(ncfile, 'lat');
    time_axis = ncread(ncfile, 'time');
end

function time_dt = netcdf_time_to_datetime(ncfile, time_axis)
    units = string(ncreadatt(ncfile, 'time', 'units'));
    units_lower = lower(units);
    if ~contains(units_lower, 'since')
        error('Unsupported time units format: %s', units);
    end

    parts = split(units, 'since');
    unit_token = strtrim(lower(parts(1)));
    ref_str = strtrim(parts(2));
    ref_time = parse_reference_time(ref_str);

    switch unit_token
        case {'days', 'day'}
            dt = days(double(time_axis(:)));
        case {'hours', 'hour', 'hrs', 'hr'}
            dt = hours(double(time_axis(:)));
        case {'minutes', 'minute', 'mins', 'min'}
            dt = minutes(double(time_axis(:)));
        case {'seconds', 'second', 'secs', 'sec'}
            dt = seconds(double(time_axis(:)));
        otherwise
            error('Unsupported time unit token in NetCDF file: %s', unit_token);
    end

    time_dt = ref_time + dt;
end

function [time_indices, time_label] = select_time_indices(time_dt, temporal_mode, inst_date, range_start, range_end)
    switch lower(string(temporal_mode))
        case "instantaneous"
            [~, idx] = min(abs(time_dt - inst_date));
            time_indices = idx;
            time_label = sprintf('instantaneous %s', datestr(time_dt(idx), 'yyyy-mm-dd HH:MM:SS UTC'));
        case "timerange"
            time_indices = find(time_dt >= range_start & time_dt <= range_end);
            if isempty(time_indices)
                error('No timestamps found inside range [%s, %s].', string(range_start), string(range_end));
            end
            time_label = sprintf('time-mean %s to %s', datestr(time_dt(time_indices(1)), 'yyyy-mm-dd'), datestr(time_dt(time_indices(end)), 'yyyy-mm-dd'));
        otherwise
            error('Unsupported temporal_mode: %s', temporal_mode);
    end
end

function [field_data, selected_percentile] = read_thresholded_field_at_percentile(ncfile, varname, target_percentile)
    percentiles = double(ncread(ncfile, 'percentile'));
    [~, pidx] = min(abs(percentiles - target_percentile));
    selected_percentile = percentiles(pidx);

    raw_values = ncread(ncfile, varname);
    dim_names = get_var_dim_names(ncfile, varname);
    values = permute_to_lon_lat_time_percentile(raw_values, dim_names);
    field_data = squeeze(values(:, :, :, pidx));
end

function dim_names = get_var_dim_names(ncfile, varname)
    info = ncinfo(ncfile, varname);
    dim_names = strings(numel(info.Dimensions), 1);
    for i = 1:numel(info.Dimensions)
        dim_names(i) = string(info.Dimensions(i).Name);
    end
end

function values = permute_to_lon_lat_time_percentile(values_in, dim_names)
    desired = ["lon", "lat", "time", "percentile"];
    perm_order = zeros(1, numel(desired));
    for i = 1:numel(desired)
        idx = find(dim_names == desired(i), 1);
        if isempty(idx)
            error('Variable is missing required dimension "%s".', desired(i));
        end
        perm_order(i) = idx;
    end
    values = permute(values_in, perm_order);
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