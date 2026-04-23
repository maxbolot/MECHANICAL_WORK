function data_smoothed = running_mean_omitnan_cols(data_in, window_days)
    % Applies centered moving-average smoothing independently to each column.
    [ntime, ncol] = size(data_in);
    half_win = floor(window_days / 2);
    data_smoothed = nan(ntime, ncol);
    for icol = 1:ncol
        column = data_in(:, icol);
        for it = 1:ntime
            start_idx = max(1, it - half_win);
            stop_idx = min(ntime, it + half_win);
            segment = column(start_idx:stop_idx);
            valid = isfinite(segment);
            if any(valid)
                data_smoothed(it, icol) = mean(segment(valid));
            end
        end
    end
end
