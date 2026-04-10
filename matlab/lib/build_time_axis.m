function [plot_time, plot_time_label] = build_time_axis(ncfile, time)
    plot_time = double(time(:));
    plot_time_label = 'Time';

    try
        units = ncreadatt(ncfile, 'time', 'units');
        units = string(units);
        units_lower = lower(units);

        if contains(units_lower, 'since')
            parts = split(units, 'since');
            unit_token = strtrim(lower(parts(1)));
            ref_str = strtrim(parts(2));

            ref_time = parse_reference_time(ref_str);

            switch unit_token
                case {'days', 'day'}
                    dt = days(double(time(:)));
                case {'hours', 'hour', 'hrs', 'hr'}
                    dt = hours(double(time(:)));
                case {'minutes', 'minute', 'mins', 'min'}
                    dt = minutes(double(time(:)));
                case {'seconds', 'second', 'secs', 'sec'}
                    dt = seconds(double(time(:)));
                otherwise
                    dt = [];
            end

            if ~isempty(dt)
                plot_time = ref_time + dt;
                plot_time_label = sprintf('Time (%s)', units);
            end
        else
            plot_time_label = sprintf('Time (%s)', units);
        end
    catch
        plot_time_label = 'Time (raw NetCDF values)';
    end
end
