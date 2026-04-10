function avg = weighted_nanmean(values, weights)
    values = double(values(:));
    weights = double(weights(:));

    valid = ~isnan(values) & ~isnan(weights) & (weights > 0);
    if ~any(valid)
        avg = NaN;
        return;
    end

    avg = sum(values(valid) .* weights(valid)) / sum(weights(valid));
end
