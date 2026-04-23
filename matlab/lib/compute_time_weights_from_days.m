function [time_weights_days, missing_steps] = compute_time_weights_from_days(time_days, tag)
    % Builds daily time weights and counts missing expected timesteps.
    if nargin < 2
        tag = 'simulation';
    end

    ntime = numel(time_days);

    if ntime == 0
        time_weights_days = [];
        missing_steps = 0;
        return;
    end

    if any(~isfinite(time_days))
        error('Invalid non-finite time values encountered while computing %s weights.', tag);
    end

    time_days = sort(time_days(:));
    if ntime > 1
        dt = diff(time_days);
        if any(dt <= 0)
            error('Time axis is not strictly increasing after sorting (%s).', tag);
        end
        if median(dt) > 30
            error('Unrealistic median timestep (%g days) in %s. Time conversion likely failed.', median(dt), tag);
        end
    end

    time_weights_days = ones(ntime, 1);

    if ntime <= 1
        missing_steps = 0;
        return;
    end

    tol_hours = get_missing_timestep_tolerance_hours_local();
    tol_days = tol_hours / 24.0;
    missing_steps = count_missing_daily_schedule_local(time_days, tol_days);
end


function tol_hours = get_missing_timestep_tolerance_hours_local()
    tol_hours = 18.0;
    tol_str = getenv('MISSING_TIMESTEP_TOL_HOURS');
    if isempty(tol_str)
        return;
    end

    tol_val = str2double(tol_str);
    if ~isfinite(tol_val) || (tol_val <= 0)
        error('MISSING_TIMESTEP_TOL_HOURS must be positive numeric. Got: %s', tol_str);
    end
    tol_hours = tol_val;
end


function missing_steps = count_missing_daily_schedule_local(time_days, tol_days)
    if numel(time_days) <= 1
        missing_steps = 0;
        return;
    end

    expected = (time_days(1):1.0:time_days(end)).';
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
    end
end
