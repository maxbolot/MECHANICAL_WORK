function labels = build_percentile_labels(percentile)
    % Creates compact percentile legend labels from threshold coordinates.
    labels = cell(numel(percentile), 1);
    for i = 1:numel(percentile)
        labels{i} = sprintf('p%g', 100.0 * double(percentile(i)));
    end
end
