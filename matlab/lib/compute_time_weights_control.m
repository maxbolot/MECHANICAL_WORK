function [time_weights_days, missing_steps] = compute_time_weights_control(time_axis, ncfile)
    % Build per-sample time weights for the control simulation schedule:
    % 5-day cadence before 2021-05-27, 1-day cadence on/after 2021-05-27.
    %
    % Outputs:
    %   time_weights_days  per-sample duration weights in days
    %   missing_steps      count of expected scheduled outputs not present
    cutoff_days = datenum(2021, 5, 27, 0, 0, 0);

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
        error('Invalid non-finite time values encountered while computing control weights.');
    end

    % Work with a sorted column vector so downstream matching is monotone.
    time_days = sort(time_days(:));
    if ntime > 1
        dt = diff(time_days);
        if any(dt <= 0)
            error('Time axis is not strictly increasing after sorting.');
        end
        % Large median spacing usually means the units were interpreted
        % incorrectly rather than the dataset genuinely having monthly output.
        if median(dt) > 30
            error('Unrealistic median timestep (%g days). Time units conversion likely failed.', median(dt));
        end
    end

    % Control schedule: 5-day segment, then 1-day segment.
    time_weights_days = ones(ntime, 1);
    % Samples before the cutoff represent 5-day means; samples on/after the
    % cutoff represent 1-day means.
    time_weights_days(time_days < cutoff_days) = 5.0;

    missing_steps = 0;
    if ntime > 1
        % Allow modest clock offsets from file metadata while still keeping
        % enough separation to distinguish neighboring scheduled outputs.
        tol_days = 3.0 / 24.0;  % 3-hour tolerance: tolerates ~1.5-hour FRE timestamp offsets while remaining << 0.5-day (half min step)

        % Prefer filename-based schedule matching from ncrcat history when
        % available. This is robust to midpoint clock-phase shifts around
        % cadence transitions.
        filename_starts_days = extract_input_start_dates_from_history(ncfile);
        if ~isempty(filename_starts_days)
            missing_steps = count_missing_from_start_dates(filename_starts_days, cutoff_days, tol_days);
            return;
        end

        % Fallback path: infer missing samples directly from timestamp axis.
        % This is less robust than filename history, but still usable when
        % the history attribute is unavailable or stripped.
        tcur = time_days(1);
        tstop = time_days(end);
        expected = tcur;
        iter_guard = 0;
        max_iter = 500000;

        % Generate the ideal schedule between the first and last observed
        % timestamps using the control cadence rules.
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

        % Two-pointer pass keeps matching linear in expected+observed length.
        % Observed times that are earlier than the current expected slot are
        % skipped until the first candidate match is reached.
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
    % Advance by the cadence active at the current sample/start time.
    if tcur_days < cutoff_days
        tnext = tcur_days + 5.0;
    else
        tnext = tcur_days + 1.0;
    end
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

    % Match only canonical single-date filenames, not repaired aliases or
    % date-range products.
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

function missing_steps = count_missing_from_start_dates(start_days, cutoff_days, tol_days)
    % Count gaps in the ideal sequence of file start times.
    if numel(start_days) <= 1
        missing_steps = 0;
        return;
    end

    % Build the ideal sequence of file start times from first to last sample.
    expected = start_days(1);
    tcur = start_days(1);
    tstop = start_days(end);
    iter_guard = 0;
    max_iter = 500000;

    while tcur < (tstop - tol_days)
        tcur = next_scheduled_time_days(tcur, cutoff_days);
        expected(end + 1, 1) = tcur; %#ok<AGROW>

        iter_guard = iter_guard + 1;
        if iter_guard > max_iter
            error('Exceeded schedule generation guard while counting missing start dates.');
        end
    end

    % Count expected starts that cannot be matched to any observed start.
    % Matching is tolerant to small clock offsets but not to full cadence gaps.
    missing_steps = 0;
    iobs = 1;
    nobs = numel(start_days);
    for iexp = 1:numel(expected)
        while (iobs <= nobs) && (start_days(iobs) < expected(iexp) - tol_days)
            iobs = iobs + 1;
        end

        if (iobs > nobs) || (abs(start_days(iobs) - expected(iexp)) > tol_days)
            missing_steps = missing_steps + 1;
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

    if all(t > 1.0e6 & t < 1.0e10)
        % Fallback for numeric YYYYMMDD.f encodings where the integer part is
        % the calendar date and the fractional part is the fraction of a day.
        day_part = floor(t);
        frac_part = t - day_part;
        time_days = zeros(size(t));
        for i = 1:numel(t)
            ymd = sprintf('%08d', round(day_part(i)));
            time_days(i) = datenum(ymd, 'yyyymmdd') + frac_part(i);
        end
        return;
    end

    error('Unable to convert time axis to day numbers for control time weighting. File: %s', ncfile);
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