function [lift3d_out, precip3d_out, common_time_days] = intersect_common_daily_time_range( ...
    % Restricts lift and precipitation fields to shared daily timestamps.
    lift3d, lift_time_axis, lift_ncfile, precip3d, precip_time_axis, precip_ncfile, tag)

    lift_time_days = normalize_time_axis_to_datenum(lift_time_axis, lift_ncfile);
    precip_time_days = normalize_time_axis_to_datenum(precip_time_axis, precip_ncfile);

    lift_day_key = floor(lift_time_days + 1.0e-9);
    precip_day_key = floor(precip_time_days + 1.0e-9);

    if numel(unique(lift_day_key)) ~= numel(lift_day_key)
        error('Lift file is not daily aggregated over unique days: %s', lift_ncfile);
    end
    if numel(unique(precip_day_key)) ~= numel(precip_day_key)
        error('Precip file is not daily aggregated over unique days: %s', precip_ncfile);
    end

    [common_day_key, lift_idx, precip_idx] = intersect(lift_day_key, precip_day_key, 'stable');
    if isempty(common_day_key)
        error('No common daily time between %s and %s.', lift_ncfile, precip_ncfile);
    end

    if (numel(common_day_key) < numel(lift_day_key)) || (numel(common_day_key) < numel(precip_day_key))
        fprintf(['Restricting %s to %d common daily samples between lift (%d) and precip (%d).\n'], ...
            tag, numel(common_day_key), numel(lift_day_key), numel(precip_day_key));
    end

    lift3d_out = lift3d(:, :, lift_idx);
    precip3d_out = precip3d(:, :, precip_idx);
    common_time_days = lift_time_days(lift_idx);
end
