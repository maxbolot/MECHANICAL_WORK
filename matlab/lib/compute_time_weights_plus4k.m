function [time_weights_days, missing_steps] = compute_time_weights_plus4k(time_axis, ncfile)
    % Build per-sample time weights for the +4K, 1270 ppmv scenario:
    % 5-day cadence before 2020-05-12, 1-day cadence on/after 2020-05-12.
    cutoff_days = datenum(2020, 5, 12, 0, 0, 0);

    if nargin < 2
        ncfile = '';
    end

    time_days = normalize_time_axis_to_datenum(time_axis, ncfile);
    ntime = numel(time_days);

    if ntime == 0
        time_weights_days = [];
        missing_steps = 0;
        return;
    end

    if any(~isfinite(time_days))
        error('Invalid non-finite time values encountered while computing weights.');
    end

    % Ensure monotonic order and reasonable spacing for schedule logic.
    time_days = sort(time_days(:));
    if ntime > 1
        dt = diff(time_days);
        if any(dt <= 0)
            error('Time axis is not strictly increasing after sorting.');
        end

        % Catch unit-conversion mistakes early instead of looping for a long time.
        if median(dt) > 30
            error('Unrealistic median timestep (%g days). Time units conversion likely failed.', median(dt));
        end
    end

    time_weights_days = ones(ntime, 1);
    time_weights_days(time_days < cutoff_days) = 5.0;

    missing_steps = 0;
    if ntime > 1
        % One-minute tolerance in days.
        tol_days = 1.0 / (24.0 * 60.0);

        % Generate the expected schedule between first and last available samples.
        tcur = time_days(1);
        tstop = time_days(end);
        expected = tcur;
        iter_guard = 0;
        max_iter = 500000;

        while tcur < (tstop - tol_days)
            tnext = next_scheduled_time_days(tcur, cutoff_days);
            if tnext <= (tcur + eps)
                error('Non-advancing expected schedule step encountered.');
            end
            expected(end + 1, 1) = tnext; %#ok<AGROW>
            tcur = tnext;

            iter_guard = iter_guard + 1;
            if iter_guard > max_iter
                error('Exceeded schedule generation guard. Check time units conversion.');
            end
        end

        % Two-pointer match between expected and observed times.
        iobs = 1;
        for iexp = 1:numel(expected)
            while (iobs <= ntime) && (time_days(iobs) < expected(iexp) - tol_days)
                iobs = iobs + 1;
            end

            if (iobs > ntime) || (abs(time_days(iobs) - expected(iexp)) > tol_days)
                missing_steps = missing_steps + 1;
            end
        end
    end
end

function tnext = next_scheduled_time_days(tcur_days, cutoff_days)
    if tcur_days < cutoff_days
        tnext = tcur_days + 5.0;
    else
        tnext = tcur_days + 1.0;
    end
end

function time_days = normalize_time_axis_to_datenum(time_axis, ncfile)
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
                ref_days = parse_reference_time_datenum(ref_str);

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

    % Handle "day as %Y%m%d.%f" convention if present as numeric.
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

    error('Unable to convert time axis to day numbers for time weighting. File: %s', ncfile);
end

function ref_days = parse_reference_time_datenum(ref_str)
    ref_str = strrep(strtrim(char(ref_str)), 'T', ' ');
    formats = {'yyyy-mm-dd HH:MM:SS', 'yyyy-mm-dd HH:MM', 'yyyy-mm-dd'};

    for i = 1:numel(formats)
        try
            ref_days = datenum(ref_str, formats{i});
            if isfinite(ref_days)
                return;
            end
        catch
        end
    end

    ref_days = datenum(ref_str);
    if ~isfinite(ref_days)
        error('Unable to parse reference time from time units: %s', ref_str);
    end
end
