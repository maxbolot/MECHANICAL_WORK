function [plot_time, plot_time_label] = build_time_axis(ncfile, time)
    % Converts raw time coordinates into plottable vectors and labels.
    plot_time = double(time(:));
    plot_time_label = 'Time';

    try
        units = ncreadatt(ncfile, 'time', 'units');
        units = string(units);
        units_lower = lower(units);

        if contains(units_lower, 'since')
            plot_time = netcdf_time_to_datetime(time, ncfile);
            plot_time_label = sprintf('Time (%s)', units);
        else
            plot_time_label = sprintf('Time (%s)', units);
        end
    catch
        plot_time_label = 'Time (raw NetCDF values)';
    end
end
