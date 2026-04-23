function values = read_thresholded_var(ncfile, varname)
    % Reads and reorders thresholded variables to [lon, lat, time, percentile].
    try
        raw_values = ncread(ncfile, varname);
        info = ncinfo(ncfile, varname);
        dim_names = strings(numel(info.Dimensions), 1);
        for i = 1:numel(info.Dimensions)
            dim_names(i) = string(info.Dimensions(i).Name);
        end
        values = permute_to_lon_lat_time_percentile(raw_values, dim_names);
    catch ME
        error('Error reading %s variable: %s\nAvailable variables: %s', varname, ME.message, get_available_vars(ncfile));
    end
end
