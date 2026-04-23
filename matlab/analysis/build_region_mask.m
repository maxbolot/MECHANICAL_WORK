function mask = build_region_mask(lat, lon, region)
    % Builds lon/lat index masks and latitude weights for configured regions.
    lat_values = double(lat(:));
    lon_values = double(lon(:));
    lon_wrapped = mod(lon_values, 360);

    if isempty(region.lat_bounds)
        lat_idx = 1:numel(lat_values);
    else
        lat_idx = find(lat_values >= region.lat_bounds(1) & lat_values <= region.lat_bounds(2));
    end

    if isempty(region.lon_bounds)
        lon_idx = 1:numel(lon_values);
    else
        lon_idx = find(lon_wrapped >= region.lon_bounds(1) & lon_wrapped <= region.lon_bounds(2));
    end

    if isempty(lat_idx) || isempty(lon_idx)
        error('No grid points found for region %s.', region.name);
    end

    lat_selected = double(lat(lat_idx));
    lat_weights = cosd(lat_selected(:));
    if all(lat_weights == 0)
        error('Latitude weights are all zero in the selected %s region.', region.name);
    end

    switch region.name
        case 'global'
            selection_message = 'Global region selected: all latitudes and longitudes';
        case 'tropics'
            selection_message = sprintf('Tropical latitude band selected: %g to %g degrees', min(lat(lat_idx)), max(lat(lat_idx)));
        case 'maritime_continent'
            selection_message = sprintf('Maritime continent region selected: lat [%g, %g], lon [%g, %g] (degrees east)', ...
                min(lat(lat_idx)), max(lat(lat_idx)), min(lon_wrapped(lon_idx)), max(lon_wrapped(lon_idx)));
        case 'northern_midlatitudes'
            selection_message = sprintf('Northern midlatitudes selected: lat [%g, %g], all longitudes', ...
                min(lat(lat_idx)), max(lat(lat_idx)));
        case 'southern_midlatitudes'
            selection_message = sprintf('Southern midlatitudes selected: lat [%g, %g], all longitudes', ...
                min(lat(lat_idx)), max(lat(lat_idx)));
        case 'southern_poles'
            selection_message = sprintf('Southern poles selected: lat [%g, %g], all longitudes', ...
                min(lat(lat_idx)), max(lat(lat_idx)));
        otherwise
            selection_message = sprintf('Region selected: %s', region.name);
    end

    mask.lat_idx = lat_idx;
    mask.lon_idx = lon_idx;
    mask.lat_weights = lat_weights;
    mask.selection_message = selection_message;
end
