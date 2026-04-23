function mean_map = weighted_mean_over_time_map(values3d, time_weights_days)
    % Computes time-weighted map means while handling NaN coverage gaps.
    weights = reshape(double(time_weights_days(:)), 1, 1, []);
    valid = ~isnan(values3d);

    weighted_num = sum(values3d .* weights, 3, 'omitnan');
    weighted_den = sum(double(valid) .* weights, 3, 'omitnan');

    mean_map = weighted_num ./ weighted_den;
    mean_map(weighted_den <= 0) = NaN;
end
