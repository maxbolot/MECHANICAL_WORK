function stats = compute_thresholded_bin_stats(cfg)
% COMPUTE_THRESHOLDED_BIN_STATS  Compute per-percentile-rank-bin averages of
%   work and lift for a given scenario/region configuration.
%
%   stats = compute_thresholded_bin_stats(cfg)
%
%   cfg must contain:
%     cfg.scenario  — from scenario_control() / scenario_plus4k()
%     cfg.region    — from region_tropics() / etc.
%
%   Returned struct fields:
%     bin_lower, bin_upper, bin_mass, bin_center  (nbin x 1)
%     work_bin_avg, lift_bin_avg, ratio_bin        (nbin x 1)
%     labels                                        (nbin x 1 cell)

    ncfile = cfg.scenario.thresholded_ncfile;
    ensure_input_exists(ncfile);

    [lat, lon, ~, ~] = read_thresholded_coordinates(ncfile);

    mask = build_region_mask(lat, lon, cfg.region);

    ts = compute_thresholded_region_spatial_timeseries(ncfile, mask);
    time = ts.time;
    percentile = ts.percentile;
    ntime = ts.ntime;
    nperc = ts.nperc;
    work_spatial_avg = ts.work_spatial_avg;
    lift_spatial_avg = ts.lift_spatial_avg;

    [time_weights_days, ~] = cfg.scenario.time_weight_fn(time, ncfile);
    if numel(time_weights_days) ~= ntime
        error('Time-weight vector length (%d) does not match data time length (%d).', numel(time_weights_days), ntime);
    end

    work_avg = nan(nperc, 1);
    lift_avg = nan(nperc, 1);
    for ip = 1:nperc
        work_avg(ip) = weighted_nanmean(work_spatial_avg(:, ip), time_weights_days);
        lift_avg(ip) = weighted_nanmean(lift_spatial_avg(:, ip), time_weights_days);
    end

    [work_total_avg, lift_total_avg] = compute_unthresholded_region_time_means(cfg, mask);

    stats = build_percentile_bin_stats(percentile, work_avg, lift_avg, work_total_avg, lift_total_avg);
end
