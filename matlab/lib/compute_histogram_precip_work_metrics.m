function out = compute_histogram_precip_work_metrics(hist_data, time_weights_days)
    if numel(time_weights_days) ~= hist_data.ntime
        error('Time-weight vector length (%d) does not match data time length (%d).', numel(time_weights_days), hist_data.ntime);
    end

    % 1) Area-weighted precipitation frequency
    hist_area_sum = sum(hist_data.hist_area, 2, 'omitnan');
    hist_area_freq = hist_data.hist_area ./ hist_area_sum;
    hist_area_freq(hist_area_sum <= 0, :) = NaN;
    out.pr_freq_timeavg = weighted_time_mean_rows(hist_area_freq, time_weights_days);

    % 2) Precipitation-binned lift/work from weighted first moments
    work_centers_bcast = reshape(hist_data.work_centers, 1, [], 1);

    lift_numerator = squeeze(sum(hist_data.hist2d_lift .* work_centers_bcast, 2, 'omitnan'));
    lift_denominator = squeeze(sum(hist_data.hist2d_lift, 2, 'omitnan'));
    work_numerator = squeeze(sum(hist_data.hist2d_work .* work_centers_bcast, 2, 'omitnan'));
    work_denominator = squeeze(sum(hist_data.hist2d_work, 2, 'omitnan'));

    w_col = time_weights_days(:);

    lift_num_w = lift_numerator .* w_col;
    lift_num_w(lift_denominator <= 0) = 0;
    lift_den_w = lift_denominator .* w_col;
    lift_den_w(lift_denominator <= 0) = 0;

    work_num_w = work_numerator .* w_col;
    work_num_w(work_denominator <= 0) = 0;
    work_den_w = work_denominator .* w_col;
    work_den_w(work_denominator <= 0) = 0;

    out.lift_by_pr_timeavg = sum(lift_num_w, 1) ./ sum(lift_den_w, 1);
    out.work_by_pr_timeavg = sum(work_num_w, 1) ./ sum(work_den_w, 1);
    out.lift_by_pr_timeavg(sum(lift_den_w, 1) <= 0) = NaN;
    out.work_by_pr_timeavg(sum(work_den_w, 1) <= 0) = NaN;

    out.ratio_by_pr_timeavg = out.lift_by_pr_timeavg ./ out.work_by_pr_timeavg;
end
