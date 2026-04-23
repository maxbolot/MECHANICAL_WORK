function out = compute_thresholded_region_spatial_timeseries(ncfile, mask)
    % Computes area-weighted thresholded regional timeseries for work and lift.
    [~, ~, time, percentile] = read_thresholded_coordinates(ncfile);

    work4d = read_thresholded_var(ncfile, 'work');
    lift4d = read_thresholded_var(ncfile, 'lift');

    ntime = size(work4d, 3);
    nperc = size(work4d, 4);
    if nperc ~= numel(percentile)
        error('Percentile dimension mismatch: data=%d, coordinate=%d', nperc, numel(percentile));
    end

    work_region = work4d(mask.lon_idx, mask.lat_idx, :, :);
    lift_region = lift4d(mask.lon_idx, mask.lat_idx, :, :);
    weight_4d = reshape(mask.lat_weights, 1, [], 1, 1);

    work_num = squeeze(sum(work_region .* weight_4d, [1, 2], 'omitnan'));
    work_den = squeeze(sum(double(~isnan(work_region)) .* weight_4d, [1, 2], 'omitnan'));
    work_spatial_avg = orient_time_percentile(work_num ./ work_den, ntime, nperc);

    lift_num = squeeze(sum(lift_region .* weight_4d, [1, 2], 'omitnan'));
    lift_den = squeeze(sum(double(~isnan(lift_region)) .* weight_4d, [1, 2], 'omitnan'));
    lift_spatial_avg = orient_time_percentile(lift_num ./ lift_den, ntime, nperc);

    out.time = time;
    out.percentile = percentile;
    out.work4d_size = size(work4d);
    out.lift4d_size = size(lift4d);
    out.ntime = ntime;
    out.nperc = nperc;
    out.work_spatial_avg = work_spatial_avg;
    out.lift_spatial_avg = lift_spatial_avg;
end
