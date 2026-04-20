function [time_weights_days, missing_steps] = compute_time_weights_plus4k(time_axis, ncfile)
    % Build per-sample time weights for the +4K, 1270 ppmv scenario.
    % Time cadence is now treated as uniformly 1 day everywhere.
    %
    % Outputs:
    %   time_weights_days  per-sample duration weights in days
    %   missing_steps      count of expected scheduled outputs not present
    if nargin < 2
        ncfile = '';
    end

    % Convert NetCDF time to MATLAB datenum so schedule logic is unit-agnostic.
    % This accepts standard "<unit> since <reference>" metadata as well as
    % the numeric YYYYMMDD.f convention seen in some workflow outputs.
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

    missing_steps = 0;
    if ntime > 1
        % Allow larger clock offsets in timestamp metadata. Override with:
        %   export MISSING_TIMESTEP_TOL_HOURS=<positive number>
        tol_hours = get_missing_timestep_tolerance_hours();
        tol_days = tol_hours / 24.0;

        % Use the NetCDF time axis as the authoritative source for cadence
        % checks. History strings from concat commands are often truncated and
        % can inflate missing-step counts.
        missing_steps = count_missing_daily_schedule(time_days, tol_days);
    end
end

function tol_hours = get_missing_timestep_tolerance_hours()
    % Default to an 18-hour tolerance to reduce sensitivity to timestamp
    % phase shifts while still being stricter than a full day.
    tol_hours = 18.0;

    tol_str = getenv('MISSING_TIMESTEP_TOL_HOURS');
    if isempty(tol_str)
        return;
    end

    tol_val = str2double(tol_str);
    if ~isfinite(tol_val) || (tol_val <= 0)
        error('MISSING_TIMESTEP_TOL_HOURS must be a positive numeric value. Got: %s', tol_str);
    end
    tol_hours = tol_val;
end

function start_days = extract_input_start_dates_from_history(ncfile)
    % Recover canonical input dates from the ncrcat command stored in the
    % global history attribute of concatenated output files.
    start_days = [];
    if isempty(ncfile)
        return;
    end

    try
        hist = char(ncreadatt(ncfile, '/', 'history'));
    catch
        return;
    end

    % Match canonical input filenames like work_YYYYMMDDHH.nc.
    tokens = regexp(hist, 'work_([0-9]{10})\.nc', 'tokens');
    if isempty(tokens)
        return;
    end

    % Keep first occurrence order to preserve original concat ordering.
    labels = strings(numel(tokens), 1);
    for i = 1:numel(tokens)
        labels(i) = string(tokens{i}{1});
    end
    labels = unique(labels, 'stable');

    % Convert YYYYMMDDHH labels to MATLAB datenums for schedule matching.
    start_days = zeros(numel(labels), 1);
    for i = 1:numel(labels)
        v = char(labels(i));
        y = str2double(v(1:4));
        m = str2double(v(5:6));
        d = str2double(v(7:8));
        h = str2double(v(9:10));
        start_days(i) = datenum(y, m, d, h, 0, 0);
    end
    start_days = sort(start_days);
end

function missing_steps = count_missing_daily_schedule(time_days, tol_days)
    % Count gaps in an ideal daily sequence spanning the observed range.
    time_days = sort(time_days(:));
    if numel(time_days) <= 1
        missing_steps = 0;
        return;
    end

    expected = time_days(1):1.0:time_days(end);
    expected = expected(:);
    iter_guard = 0;
    max_iter = 500000;

    missing_steps = 0;
    iobs = 1;
    nobs = numel(time_days);
    for iexp = 1:numel(expected)
        while (iobs <= nobs) && (time_days(iobs) < expected(iexp) - tol_days)
            iobs = iobs + 1;
        end

        if (iobs > nobs) || (abs(time_days(iobs) - expected(iexp)) > tol_days)
            missing_steps = missing_steps + 1;
        end

        iter_guard = iter_guard + 1;
        if iter_guard > max_iter
            error('Exceeded daily schedule generation guard while counting missing times.');
        end
    end
end

function time_days = normalize_time_axis_to_datenum(time_axis, ncfile)
    % Normalize supported time encodings into MATLAB datenum values.
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
                % Parse CF-style units such as "days since ..." or
                % "minutes since ...".
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
        % The integer part is the date; the fractional part is the fraction
        % of a day after 00:00 UTC.
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
    % Accept a few common timestamp formats before falling back to datenum's
    % generic parser.
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
