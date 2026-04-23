function [lat, lon, time, percentile] = read_thresholded_coordinates(ncfile)
    % Reads lat/lon/time/percentile coordinates from a thresholded NetCDF file.
    try
        lat = ncread(ncfile, 'lat');
        lon = ncread(ncfile, 'lon');
        time = ncread(ncfile, 'time');
        percentile = ncread(ncfile, 'percentile');
    catch ME
        error('Error reading thresholded coordinates: %s', ME.message);
    end
end
