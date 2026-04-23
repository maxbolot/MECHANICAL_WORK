function time_days = normalize_time_axis_to_datenum(time_axis, ncfile)
    % Normalizes supported time encodings into MATLAB datenum day numbers.
    if nargin < 2
        ncfile = '';
    end

    if isdatetime(time_axis)
        time_days = datenum(time_axis(:));
        return;
    end

    t = double(time_axis(:));

    if ~isempty(ncfile)
        try
            units = string(ncreadatt(ncfile, 'time', 'units'));
            units_lower = lower(units);

            if contains(units_lower, 'since')
                parts = split(units, 'since');
                unit_token = strtrim(lower(parts(1)));
                ref_str = strtrim(parts(2));
                ref_days = datenum(parse_reference_time(ref_str));

                switch unit_token
                    case {'days', 'day'}
                        time_days = ref_days + t;
                        return;
                    case {'hours', 'hour', 'hrs', 'hr'}
                        time_days = ref_days + t / 24.0;
                        return;
                    case {'minutes', 'minute', 'mins', 'min'}
                        time_days = ref_days + t / (24.0 * 60.0);
                        return;
                    case {'seconds', 'second', 'secs', 'sec'}
                        time_days = ref_days + t / (24.0 * 60.0 * 60.0);
                        return;
                end
            end
        catch
        end
    end

    if all(t > 1.0e6 & t < 1.0e10)
        day_part = floor(t);
        frac_part = t - day_part;
        time_days = zeros(size(t));
        for i = 1:numel(t)
            ymd = sprintf('%08d', round(day_part(i)));
            time_days(i) = datenum(ymd, 'yyyymmdd') + frac_part(i);
        end
        return;
    end

    error('Unable to convert time axis to day numbers. File: %s', ncfile);
end
