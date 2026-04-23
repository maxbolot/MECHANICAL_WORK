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

    time_days = normalize_time_axis_to_datenum(time_axis, ncfile);
    [time_weights_days, missing_steps] = compute_time_weights_from_days(time_days, 'plus4k');
end
