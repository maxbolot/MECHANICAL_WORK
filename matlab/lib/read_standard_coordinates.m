function [lat, lon, time] = read_standard_coordinates(ncfile)
    % Reads canonical lat/lon/time coordinates from a standard NetCDF file.
    try
        lat = ncread(ncfile, 'lat');
        lon = ncread(ncfile, 'lon');
        time = ncread(ncfile, 'time');
    catch ME
        error('Error reading standard coordinates: %s', ME.message);
    end
end
