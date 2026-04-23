function report = regression_check_thresholded_bin_stats(cfg)
% REGRESSION_CHECK_THRESHOLDED_BIN_STATS Compare thresholded bin stats from
% run_thresholded_region_analysis() and compute_thresholded_bin_stats().

    if nargin < 1 || isempty(cfg)
        cfg = struct();
        cfg.scenario = scenario_control();
        cfg.region = region_tropics();
    end

    if ~isfield(cfg, 'scenario')
        cfg.scenario = scenario_control();
    end
    if ~isfield(cfg, 'region')
        cfg.region = region_tropics();
    end

    cfg_for_run = cfg;
    cfg_for_run.smooth_before_plotting = false;
    cfg_for_run.smooth_window_days = 5;
    cfg_for_run.make_plots = false;
    cfg_for_run.print_summary = false;

    run_out = run_thresholded_region_analysis(cfg_for_run);
    stats_direct = compute_thresholded_bin_stats(cfg);

    report = struct();
    report.scenario_name = cfg.scenario.name;
    report.region_name = cfg.region.result_label;
    report.max_abs_diff_work_bin = max(abs(run_out.bin_stats.work_bin_avg - stats_direct.work_bin_avg), [], 'omitnan');
    report.max_abs_diff_lift_bin = max(abs(run_out.bin_stats.lift_bin_avg - stats_direct.lift_bin_avg), [], 'omitnan');
    report.max_abs_diff_ratio_bin = max(abs(run_out.bin_stats.ratio_bin - stats_direct.ratio_bin), [], 'omitnan');

    if isempty(report.max_abs_diff_work_bin)
        report.max_abs_diff_work_bin = 0;
    end
    if isempty(report.max_abs_diff_lift_bin)
        report.max_abs_diff_lift_bin = 0;
    end
    if isempty(report.max_abs_diff_ratio_bin)
        report.max_abs_diff_ratio_bin = 0;
    end

    fprintf('\nRegression check (thresholded bin stats)\n');
    fprintf('  scenario: %s\n', report.scenario_name);
    fprintf('  region  : %s\n', report.region_name);
    fprintf('  max |work_bin_avg diff| : %.12g\n', report.max_abs_diff_work_bin);
    fprintf('  max |lift_bin_avg diff| : %.12g\n', report.max_abs_diff_lift_bin);
    fprintf('  max |ratio_bin diff|    : %.12g\n', report.max_abs_diff_ratio_bin);
end
