function [var_out, dim_names] = read_var_as_lon_lat_time(ncfile, varname)
    % Reads and reorders a variable to [lon, lat, time] based on named dimensions.
    info = ncinfo(ncfile, varname);
    dim_names = lower(string({info.Dimensions.Name}));

    raw = double(ncread(ncfile, varname));

    lon_pos = find(dim_names == "lon", 1);
    lat_pos = find(dim_names == "lat", 1);
    time_pos = find(dim_names == "time", 1);

    if isempty(lon_pos) || isempty(lat_pos) || isempty(time_pos)
        error('Variable %s in %s does not contain lon/lat/time dims.', varname, ncfile);
    end

    perm = [lon_pos, lat_pos, time_pos];
    var_out = permute(raw, perm);
end
