function [work_avg, lift_avg] = compute_unthresholded_region_time_means(cfg, mask)
    % Computes weighted regional mean work and lift from standard files.
    ncfile = cfg.scenario.standard_ncfile;
    ensure_input_exists(ncfile);

    [~, ~, time] = read_standard_coordinates(ncfile);
    work3d = read_required_var(ncfile, 'work');
    lift3d = read_required_var(ncfile, 'lift');

    work_region = work3d(mask.lon_idx, mask.lat_idx, :);
    lift_region = lift3d(mask.lon_idx, mask.lat_idx, :);
    weight_3d = reshape(mask.lat_weights, 1, [], 1);

    work_num = sum(work_region .* weight_3d, [1, 2], 'omitnan');
    work_den = sum(double(~isnan(work_region)) .* weight_3d, [1, 2], 'omitnan');
    work_spatial_avg = squeeze(work_num ./ work_den);

    lift_num = sum(lift_region .* weight_3d, [1, 2], 'omitnan');
    lift_den = sum(double(~isnan(lift_region)) .* weight_3d, [1, 2], 'omitnan');
    lift_spatial_avg = squeeze(lift_num ./ lift_den);

    time_weights_days = cfg.scenario.time_weight_fn(time, ncfile);
    if iscell(time_weights_days)
        time_weights_days = time_weights_days{1};
    end

    work_avg = weighted_nanmean(work_spatial_avg, time_weights_days);
    lift_avg = weighted_nanmean(lift_spatial_avg, time_weights_days);
end
