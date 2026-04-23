function [var3d, selected_percentile] = read_lon_lat_time_from_var(ncfile, varname, target_percentile)
    % Reads a variable as [lon, lat, time], selecting percentile slices when present.
    if nargin < 3
        target_percentile = NaN;
    end

    info = ncinfo(ncfile, varname);
    dim_names = lower(string({info.Dimensions.Name}));

    if any(dim_names == "percentile")
        values = double(ncread(ncfile, varname));
        percentiles = double(ncread(ncfile, 'percentile'));

        if ~isfinite(target_percentile)
            pidx = 1;
        else
            [~, pidx] = min(abs(percentiles - target_percentile));
        end
        selected_percentile = percentiles(pidx);

        desired = ["lon", "lat", "time", "percentile"];
        perm_order = zeros(1, numel(desired));
        for i = 1:numel(desired)
            pos = find(dim_names == desired(i), 1);
            if isempty(pos)
                error('Variable %s in %s is missing required dimension "%s".', varname, ncfile, desired(i));
            end
            perm_order(i) = pos;
        end

        values = permute(values, perm_order);
        var3d = squeeze(values(:, :, :, pidx));
        return;
    end

    selected_percentile = NaN;
    var3d = read_var_as_lon_lat_time(ncfile, varname);
end
