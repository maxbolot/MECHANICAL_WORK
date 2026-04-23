function time_dt = netcdf_time_to_datetime(time_axis, ncfile)
    % Converts NetCDF time coordinates to timezone-aware UTC datetimes.
    time_days = normalize_time_axis_to_datenum(time_axis, ncfile);
    time_dt = datetime(time_days, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
end
