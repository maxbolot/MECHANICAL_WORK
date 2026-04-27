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
        build_main_title(cfg.run_mode, control_maps, warming_maps), cfg);
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
    if ~isfield(cfg, 'wet_day_threshold')
        cfg.wet_day_threshold = 1.0e-5;
    end
    if ~isfield(cfg, 'g')
        cfg.g = 9.81;
    end
    if ~isfield(cfg, 'clim_hp')
        cfg.clim_hp = [-500, 2000];
    end
    if ~isfield(cfg, 'clim_precip')
        cfg.clim_precip = [-1.0e-4, 1.0e-4];
    end
    if ~isfield(cfg, 'contour_levels_hp')
        cfg.contour_levels_hp = [3500, 4000, 5000];
    end
    if ~isfield(cfg, 'contour_levels_precip')
        cfg.contour_levels_precip = [2, 5, 10, 15, 20] .* 1.0e-6;
    end
    if ~isfield(cfg, 'show_contour_labels')
        cfg.show_contour_labels = false;
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

    [lift3d, precip3d, common_time_days, lift_idx] = intersect_common_daily_time_range( ...
        lift3d, lift_time_axis, files.lift_ncfile, ...
        precip3d, precip_time_axis, files.precip_ncfile, ...
        tag);

    if is_wet_day_mode_local(cfg.run_mode)
        wet_mask = (precip3d > cfg.wet_day_threshold);
        lift3d(~wet_mask) = NaN;
        precip3d(~wet_mask) = NaN;
    end

    [time_weights_days, missing_steps] = compute_time_weights_from_days(common_time_days, tag);

    out.lon = lon;
    out.lat = lat;
    if is_thresholded_mode_local(cfg.run_mode)
        [event_count3d, ~] = read_lon_lat_time_from_var(files.lift_ncfile, 'event_count', cfg.threshold_percentile);
        event_count3d = event_count3d(:, :, lift_idx);
        event_weights = event_count3d .* reshape(double(time_weights_days(:)), 1, 1, []);
        out.lift_mean = weighted_mean_over_time_map_with_cell_weights_local(lift3d, event_weights);
        out.precip_mean = weighted_mean_over_time_map_with_cell_weights_local(precip3d, event_weights);
    else
        out.lift_mean = weighted_mean_over_time_map(lift3d, time_weights_days);
        out.precip_mean = weighted_mean_over_time_map(precip3d, time_weights_days);
    end
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
        case "wet_day_only"
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
        case "prate_thresholded_by_lat_band"
            if isfield(scenario, 'thresholded_by_lat_band_ncfile') && isfield(scenario, 'thresholded_by_lat_band_precip_ncfile')
                files.lift_ncfile = scenario.thresholded_by_lat_band_ncfile;
                files.precip_ncfile = scenario.thresholded_by_lat_band_precip_ncfile;
            else
                error(['Scenario preset does not define explicit lat-band thresholded files. ', ...
                    'Expected fields: thresholded_by_lat_band_ncfile and thresholded_by_lat_band_precip_ncfile.']);
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


function tf = is_thresholded_mode_local(run_mode)
    rm = lower(string(run_mode));
    tf = (rm == "prate_thresholded") || (rm == "prate_thresholded_by_lat_band");
end


function tf = is_wet_day_mode_local(run_mode)
    rm = lower(string(run_mode));
    tf = (rm == "wet_day_only");
end


function mean_map = weighted_mean_over_time_map_with_cell_weights_local(values3d, weights3d)
    valid = isfinite(values3d) & isfinite(weights3d) & (weights3d > 0);
    weighted_num = sum(values3d .* weights3d .* double(valid), 3, 'omitnan');
    weighted_den = sum(weights3d .* double(valid), 3, 'omitnan');

    mean_map = weighted_num ./ weighted_den;
    mean_map(weighted_den <= 0) = NaN;
end


function plot_delta_maps(lon_in, lat_in, delta_hp, hp_control, delta_p, p_control, main_title, cfg)
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
    plot_delta_panel(ax1, X, Y, delta_hp_plot, hp_control_plot, cm, ...
        'a. \DeltaH_p with control H_p contours', 'm', cfg.clim_hp, cfg.contour_levels_hp, true, cfg.show_contour_labels);

    ax2 = nexttile(2);
    plot_delta_panel(ax2, X, Y, delta_p_plot, p_control_plot, cm, ...
        'b. \DeltaP with control P contours', 'kg m^{-2} s^{-1}', cfg.clim_precip, cfg.contour_levels_precip, false, cfg.show_contour_labels);

    sgtitle(t, main_title, 'Interpreter', 'none');
end


function plot_delta_panel(ax, X, Y, delta_field, control_field, cm, panel_title, colorbar_label, clim_override, contour_levels_override, smooth_control_contours, show_contour_labels)
    axes(ax);
    nan_background = [0.7, 0.7, 0.7];
    panel_clim = resolve_panel_clim(delta_field, clim_override);
    rgb_field = map_field_to_rgb_with_nan_local(delta_field, cm, panel_clim, nan_background);
    contour_field = control_field;
    if smooth_control_contours
        contour_field = imgaussfilt(control_field, 1.5);
    end

    axesm('mercator', 'Frame', 'on', 'Grid', 'on', 'MapLatLimit', [-60 60], ...
        'MapLonLimit', [-180 180], 'Origin', [0 0 0], ...
        'meridianlabel', 'on', 'parallellabel', 'on', 'labelformat', 'signed', ...
        'mlabellocation', 60, 'plinelocation', 20, 'mlinelocation', 45, ...
        'MLabelParallel', 'south', 'FontSize', 8);

    geoshow(Y, X, rgb_field, 'DisplayType', 'texturemap');
    shading flat;

    h_coast = gobjects(1);
    try
        h_coast = geoshow('landareas.shp', 'FaceColor', 'none', 'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
    catch
        warning('landareas.shp not found. Continuing without coastlines.');
    end

    hold on;
    contour_levels = resolve_contour_levels(contour_field, contour_levels_override);
    print_contour_levels(panel_title, contour_levels, contour_levels_override);
    if ~isempty(contour_levels)
        [cmat, ~] = contourm(Y, X, contour_field, contour_levels, 'Color', [0.2, 0.2, 0.2], 'LineWidth', 0.8);
        if show_contour_labels
            if has_labelable_contour_segments_local(cmat)
                try
                    clabel(cmat, 'Color', [0.2, 0.2, 0.2], 'FontSize', 8);
                catch ME
                    warning('Contour labeling skipped for panel "%s": %s', panel_title, ME.message);
                end
            else
                fprintf('Contour labeling skipped for panel "%s": no labelable contour segments.\n', panel_title);
            end
        end
    end
    hold off;

    axis off;
    set(gca, 'Color', 'w');
    set(gca, 'Layer', 'top');
    colormap(gca, cm);

    clim(panel_clim);

    title(panel_title, 'FontWeight', 'normal');
    ax.TitleHorizontalAlignment = 'left';

    cb = colorbar('eastoutside');
    cb.Label.String = colorbar_label;
end


function rgb_field = map_field_to_rgb_with_nan_local(field, cm, clim_vals, nan_color)
    ncolors = size(cm, 1);
    rgb_field = zeros([size(field), 3]);

    finite_mask = isfinite(field);
    rgb_field(:, :, 1) = nan_color(1);
    rgb_field(:, :, 2) = nan_color(2);
    rgb_field(:, :, 3) = nan_color(3);

    if ~any(finite_mask(:))
        return;
    end

    scaled = (field - clim_vals(1)) ./ (clim_vals(2) - clim_vals(1));
    scaled = min(max(scaled, 0), 1);
    color_idx = floor(scaled * (ncolors - 1)) + 1;

    for ichan = 1:3
        channel = rgb_field(:, :, ichan);
        channel(finite_mask) = cm(color_idx(finite_mask), ichan);
        rgb_field(:, :, ichan) = channel;
    end
end


function print_contour_levels(panel_title, contour_levels, contour_levels_override)
    if isempty(contour_levels)
        fprintf('%s contour levels: none\n', panel_title);
        return;
    end

    if isempty(contour_levels_override)
        source = 'auto';
    else
        source = 'manual';
    end

    fprintf('%s contour levels (%s): %s\n', panel_title, source, mat2str(contour_levels, 6));
end


function contour_levels = resolve_contour_levels(control_field, contour_levels_override)
    if ~isempty(contour_levels_override)
        if ~(isnumeric(contour_levels_override) && isvector(contour_levels_override) && all(isfinite(contour_levels_override(:))))
            error('contour level override must be a finite numeric vector.');
        end
        contour_levels = unique(sort(double(contour_levels_override(:).')));
        if isempty(contour_levels)
            error('contour level override cannot be empty after filtering.');
        end
        return;
    end

    contour_levels = contour_levels_for_field(control_field, 9);
end


function tf = has_labelable_contour_segments_local(cmat)
    if isempty(cmat) || size(cmat, 1) ~= 2 || size(cmat, 2) < 2
        tf = false;
        return;
    end

    tf = false;
    idx = 1;
    ncol = size(cmat, 2);
    while idx <= ncol
        npts = cmat(2, idx);
        if ~isfinite(npts) || (npts < 2)
            break;
        end
        if npts >= 3
            tf = true;
            return;
        end
        idx = idx + npts + 1;
    end
end


function panel_clim = resolve_panel_clim(delta_field, clim_override)
    if ~isempty(clim_override)
        if ~(isnumeric(clim_override) && numel(clim_override) == 2 && all(isfinite(clim_override(:))))
            error('clim override must be a finite 2-element numeric vector [cmin cmax].');
        end
        panel_clim = double(clim_override(:)).';
        if panel_clim(1) >= panel_clim(2)
            error('clim override must satisfy cmin < cmax.');
        end
        return;
    end

    cmax = max(abs(delta_field(:)), [], 'omitnan');
    if ~isfinite(cmax) || (cmax <= 0)
        cmax = 1;
    end
    panel_clim = [-cmax, cmax];
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

