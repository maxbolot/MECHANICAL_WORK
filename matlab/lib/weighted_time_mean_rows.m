function row_timeavg = weighted_time_mean_rows(data_time_by_bin, time_weights)
    % Compute weighted mean over time independently for each bin/column.
    % Inputs:
    %   data_time_by_bin: (time, nbins)
    %   time_weights:     (time, 1)

    data_time_by_bin = double(data_time_by_bin);
    time_weights = double(time_weights(:));

    if size(data_time_by_bin, 1) ~= numel(time_weights)
        error('Size mismatch: data has %d timesteps but weights has %d.', size(data_time_by_bin, 1), numel(time_weights));
    end

    valid = ~isnan(data_time_by_bin) & isfinite(data_time_by_bin) & (time_weights > 0);
    weighted_data = data_time_by_bin .* time_weights;
    weighted_data(~valid) = 0;

    sum_w = sum(time_weights .* valid, 1);
    sum_xw = sum(weighted_data, 1);

    row_timeavg = sum_xw ./ sum_w;
    row_timeavg(sum_w <= 0) = NaN;
end
