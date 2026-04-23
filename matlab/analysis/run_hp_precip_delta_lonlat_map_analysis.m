function run_hp_precip_delta_lonlat_map_analysis(cfg)
    % Computes control-vs-warming Hp and precipitation deltas and renders map panels.
    cfg = apply_defaults(cfg);

    control_scenario = scenario_control();
    warming_scenario = scenario_plus4k();

    control_maps = process_simulation_pair(control_scenario, cfg, 'control');
    warming_maps = process_simulation_pair(warming_scenario, cfg, 'warming');

    if numel(control_maps.lon) ~= numel(warming_maps.lon) || any(abs(control_maps.lon - warming_maps.lon) > 1.0e-10)
        error('Control and warming longitudes do not match.');
    end
    if numel(control_maps.lat) ~= numel(warming_maps.lat) || any(abs(control_maps.lat - warming_maps.lat) > 1.0e-10)
        error('Control and warming latitudes do not match.');
    end

    delta_hp = warming_maps.hp_mean - control_maps.hp_mean;
    delta_precip = warming_maps.precip_mean - control_maps.precip_mean;

    plot_delta_maps(control_maps.lon, control_maps.lat, ...
        delta_hp, control_maps.hp_mean, ...
        delta_precip, control_maps.precip_mean, ...
        build_main_title(cfg.run_mode, control_maps, warming_maps));
end


function cfg = apply_defaults(cfg)
    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    if ~isfield(cfg, 'run_mode')
        cfg.run_mode = 'non_thresholded';
    end
    if ~isfield(cfg, 'threshold_percentile')
        cfg.threshold_percentile = 0.5;
    end
    if ~isfield(cfg, 'g')
        cfg.g = 9.81;
    end
end


function out = process_simulation_pair(scenario, cfg, tag)
    files = resolve_files_for_mode(scenario, cfg.run_mode);

    ensure_file_exists(files.lift_ncfile, sprintf('%s lift/work file', tag));
    if isempty(files.precip_ncfile) || ~isfile(files.precip_ncfile)
        fprintf(['%s precip file is missing in preset for mode "%s". ' ...
            'Falling back to lift/work file for precipitation variable lookup.\n'], tag, cfg.run_mode);
        files.precip_ncfile = files.lift_ncfile;
    end

    lon = double(ncread(files.lift_ncfile, 'lon'));
    lat = double(ncread(files.lift_ncfile, 'lat'));
    lift_time_axis = double(ncread(files.lift_ncfile, 'time'));
    precip_time_axis = double(ncread(files.precip_ncfile, 'time'));

    [lift3d, lift_percentile] = read_lon_lat_time_from_var(files.lift_ncfile, 'lift', cfg.threshold_percentile);
    precip_varname = find_precip_varname(files.precip_ncfile);
    [precip3d, precip_percentile] = read_lon_lat_time_from_var(files.precip_ncfile, precip_varname, cfg.threshold_percentile);

    [lift3d, precip3d, common_time_days] = intersect_common_daily_time_range( ...
        lift3d, lift_time_axis, files.lift_ncfile, ...
        precip3d, precip_time_axis, files.precip_ncfile, ...
        tag);

    [time_weights_days, missing_steps] = compute_time_weights_from_days(common_time_days, tag);

    out.lon = lon;
    out.lat = lat;
    out.lift_mean = weighted_mean_over_time_map(lift3d, time_weights_days);
    out.precip_mean = weighted_mean_over_time_map(precip3d, time_weights_days);
    out.hp_mean = out.lift_mean ./ (cfg.g .* out.precip_mean);
    out.hp_mean(abs(out.precip_mean) < eps) = NaN;

    out.sample_count = numel(common_time_days);
    out.missing_steps = missing_steps;
    out.lift_percentile = lift_percentile;
    out.precip_percentile = precip_percentile;

    fprintf('%s mode=%s samples=%d missing=%d\n', tag, cfg.run_mode, out.sample_count, out.missing_steps);
end


function files = resolve_files_for_mode(scenario, run_mode)
    switch lower(string(run_mode))
        case "non_thresholded"
            files.lift_ncfile = scenario.standard_ncfile;
            if isfield(scenario, 'standard_precip_ncfile')
                files.precip_ncfile = scenario.standard_precip_ncfile;
            else
                files.precip_ncfile = '';
            end
        case "prate_thresholded"
            files.lift_ncfile = scenario.thresholded_ncfile;
            if isfield(scenario, 'thresholded_precip_ncfile')
                files.precip_ncfile = scenario.thresholded_precip_ncfile;
            else
                files.precip_ncfile = '';
            end
        otherwise
            error('Unsupported run_mode: %s', run_mode);
    end
end


function ensure_file_exists(ncfile, label)
    if ~isfile(ncfile)
        error('%s not found: %s', label, ncfile);
    end
end


function plot_delta_maps(lon_in, lat_in, delta_hp, hp_control, delta_p, p_control, main_title)
    [lon, delta_hp, hp_control] = normalize_and_sort_longitude(lon_in, delta_hp, hp_control);
    [~, delta_p, p_control] = normalize_and_sort_longitude(lon_in, delta_p, p_control);
    lat = double(lat_in(:));

    [lon_plot, delta_hp_plot] = append_cyclic_longitude(lon, delta_hp);
    [~, hp_control_plot] = append_cyclic_longitude(lon, hp_control);
    [~, delta_p_plot] = append_cyclic_longitude(lon, delta_p);
    [~, p_control_plot] = append_cyclic_longitude(lon, p_control);

    [X, Y] = ndgrid(lon_plot, lat);

    fig = figure('Units', 'inches', 'Position', [0, 0, 8.5, 6.5]);
    set(fig, 'Color', 'w');
    set(fig, 'DefaultAxesFontName', 'Arial');
    set(fig, 'DefaultAxesFontSize', 9);
    set(fig, 'WindowStyle', 'docked');

    t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    cm = blue_white_red_colormap(256);

    ax1 = nexttile(1);
    plot_delta_panel(ax1, X, Y, delta_hp_plot, hp_control_plot, cm, 'a. \DeltaH_p with control H_p contours', 'm');

    ax2 = nexttile(2);
    plot_delta_panel(ax2, X, Y, delta_p_plot, p_control_plot, cm, 'b. \DeltaP with control P contours', 'kg m^{-2} s^{-1}');

    sgtitle(t, main_title, 'Interpreter', 'none');
end


function plot_delta_panel(ax, X, Y, delta_field, control_field, cm, panel_title, colorbar_label)
    axes(ax);

    axesm('mercator', 'Frame', 'on', 'Grid', 'on', 'MapLatLimit', [-60 60], ...
        'MapLonLimit', [-180 180], 'Origin', [0 0 0], ...
        'meridianlabel', 'on', 'parallellabel', 'on', 'labelformat', 'signed', ...
        'mlabellocation', 60, 'plinelocation', 20, 'mlinelocation', 45, ...
        'MLabelParallel', 'south', 'FontSize', 8);

    geoshow(Y, X, delta_field, 'DisplayType', 'texturemap');
    shading flat;

    hold on;
    contour_levels = contour_levels_for_field(control_field, 9);
    if ~isempty(contour_levels)
        contourm(Y, X, control_field, contour_levels, 'k', 'LineWidth', 0.6);
    end
    hold off;

    try
        geoshow('landareas.shp', 'FaceColor', 'none', 'EdgeColor', [0.2, 0.2, 0.2], 'LineWidth', 0.5);
    catch
        warning('landareas.shp not found. Continuing without coastlines.');
    end

    axis off;
    set(gca, 'Color', 'none');
    set(gca, 'Layer', 'top');
    colormap(gca, cm);

    cmax = max(abs(delta_field(:)), [], 'omitnan');
    if ~isfinite(cmax) || (cmax <= 0)
        cmax = 1;
    end
    clim([-cmax, cmax]);

    title(panel_title, 'FontWeight', 'normal');
    ax.TitleHorizontalAlignment = 'left';

    cb = colorbar('eastoutside');
    cb.Label.String = colorbar_label;
end


function levels = contour_levels_for_field(field, nlevels)
    fmin = min(field(:), [], 'omitnan');
    fmax = max(field(:), [], 'omitnan');

    if ~isfinite(fmin) || ~isfinite(fmax) || (fmin == fmax)
        levels = [];
        return;
    end

    levels = linspace(fmin, fmax, nlevels);
end


function [lon_sorted, map1_sorted, map2_sorted] = normalize_and_sort_longitude(lon_in, map1, map2)
    lon = double(lon_in(:));
    lon = mod(lon + 180, 360) - 180;
    [lon_sorted, idx] = sort(lon);

    map1_sorted = map1(idx, :);
    map2_sorted = map2(idx, :);
end


function [lon_plot, field_plot] = append_cyclic_longitude(lon, field)
    lon_plot = [lon; lon(1) + 360.0];
    field_plot = [field; field(1, :)];
end


function cm = blue_white_red_colormap(n)
    if nargin < 1
        n = 256;
    end

    n1 = floor(n / 2);
    n2 = n - n1;

    blue = [0.1, 0.25, 0.85];
    white = [1.0, 1.0, 1.0];
    red = [0.85, 0.15, 0.1];

    x1 = linspace(0, 1, n1).';
    x2 = linspace(0, 1, n2).';

    c1 = blue + (white - blue) .* x1;
    c2 = white + (red - white) .* x2;
    cm = [c1; c2];
end


function title_text = build_main_title(run_mode, control_maps, warming_maps)
    title_text = sprintf(['Delta maps (%s): warming - control, full available 2-year period; ' ...
        'control contours overlaid (Nctrl=%d, Nwarm=%d)'], ...
        run_mode, control_maps.sample_count, warming_maps.sample_count);
end

