function summarize_work_lift(summary)
    fprintf('\n');
    fprintf('=== RESULTS ===\n');
    fprintf('Region: %s\n', summary.region_label);
    fprintf('Spatial weighting: cos(latitude)\n');
    fprintf('Time weighting: uniform 1-day cadence\n');
    fprintf('Detected missing timesteps: %d\n', summary.missing_steps);
    fprintf('Number of time steps: %d\n', summary.ntime);

    if isfield(summary, 'is_thresholded') && summary.is_thresholded
        fprintf('Number of percentile thresholds: %d\n', numel(summary.percentile_labels));
        if isfield(summary, 'smooth_before_plotting')
            if summary.smooth_before_plotting
                fprintf('Plot smoothing: enabled (%d-day running mean)\n', summary.smooth_window_days);
            else
                fprintf('Plot smoothing: disabled\n');
            end
        end
        fprintf('\n%-12s %14s %14s %14s\n', 'Percentile', 'Work', 'Lift', 'Lift/Work');
        for i = 1:numel(summary.percentile_labels)
            fprintf('%-12s %14.6f %14.6f %14.6f\n', summary.percentile_labels{i}, summary.work_avg(i), summary.lift_avg(i), summary.ratio_avg(i));
        end
        return;
    end

    fprintf('\nMechanical Work (averaged over %s and time): %.6f\n', summary.average_label, summary.work_avg);
    fprintf('Lift Work (averaged over %s and time): %.6f\n', summary.average_label, summary.lift_avg);
    fprintf('Lift/Work ratio (time-mean values): %.6f\n', summary.ratio_avg);
end