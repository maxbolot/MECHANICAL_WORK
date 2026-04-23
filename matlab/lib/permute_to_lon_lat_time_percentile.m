function values = permute_to_lon_lat_time_percentile(values_in, dim_names)
    % Reorders arrays to [lon, lat, time, percentile] dimension order.
    desired = ["lon", "lat", "time", "percentile"];
    perm_order = zeros(1, numel(desired));
    for i = 1:numel(desired)
        idx = find(dim_names == desired(i), 1);
        if isempty(idx)
            error('Variable is missing required dimension "%s". Found: %s', desired(i), strjoin(cellstr(dim_names), ', '));
        end
        perm_order(i) = idx;
    end
    values = permute(values_in, perm_order);
end
