function avg = weighted_nanmean(values, weights)
    % Computes a weighted mean while excluding non-finite samples.
    values = double(values(:));
    weights = double(weights(:));

    valid = ~isnan(values) & ~isnan(weights) & (weights > 0);
    if ~any(valid)
        avg = NaN;
        return;
    end

    avg = sum(values(valid) .* weights(valid)) / sum(weights(valid));
end
