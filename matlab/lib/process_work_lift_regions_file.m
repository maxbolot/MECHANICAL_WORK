function out = process_work_lift_regions_file(ncfile, regions, simulation_name)
    % Processes one work/lift file into region-aggregated summary diagnostics.
    lat = double(ncread(ncfile, 'lat'));
    lon = double(ncread(ncfile, 'lon'));
    time = double(ncread(ncfile, 'time'));

    [work3d, dim_names] = read_var_as_lon_lat_time(ncfile, 'work');
    [lift3d, ~] = read_var_as_lon_lat_time(ncfile, 'lift');

    fprintf('  %s work dims interpreted as [lon, lat, time] from source dims [%s]\n', ...
        ncfile, strjoin(dim_names, ', '));

    if strcmpi(simulation_name, 'warming') || strcmpi(simulation_name, 'plus4k')
        [time_weights_days, missing_steps] = compute_time_weights_plus4k(time, ncfile);
    elseif strcmpi(simulation_name, 'control')
        [time_weights_days, missing_steps] = compute_time_weights_control(time, ncfile);
    else
        error('Unknown simulation name: %s', simulation_name);
    end

    nreg = numel(regions);
    out.work = zeros(nreg, 1);
    out.lift = zeros(nreg, 1);
    out.ke = zeros(nreg, 1);
    out.ratio = zeros(nreg, 1);
    out.missing_steps = zeros(nreg, 1);

    lon_wrapped = mod(lon(:), 360);

    for ir = 1:nreg
        r = regions{ir};

        lat_idx = find(lat >= r.latmin & lat <= r.latmax);
        if isnan(r.lonmin)
            lon_idx = 1:numel(lon);
        else
            lon_idx = find(lon_wrapped >= r.lonmin & lon_wrapped <= r.lonmax);
        end

        if isempty(lat_idx) || isempty(lon_idx)
            error('No grid cells found for region: %s', r.name);
        end

        lat_weights = cosd(lat(lat_idx));
        w3 = reshape(lat_weights(:), 1, [], 1);

        wr = work3d(lon_idx, lat_idx, :);
        lr = lift3d(lon_idx, lat_idx, :);

        w_num = sum(wr .* w3, [1, 2], 'omitnan');
        w_den = sum(double(~isnan(wr)) .* w3, [1, 2], 'omitnan');
        work_ts = squeeze(w_num ./ w_den);

        l_num = sum(lr .* w3, [1, 2], 'omitnan');
        l_den = sum(double(~isnan(lr)) .* w3, [1, 2], 'omitnan');
        lift_ts = squeeze(l_num ./ l_den);

        work_avg = weighted_nanmean(work_ts, time_weights_days);
        lift_avg = weighted_nanmean(lift_ts, time_weights_days);

        out.work(ir) = work_avg;
        out.lift(ir) = lift_avg;
        out.ke(ir) = work_avg - lift_avg;
        if abs(work_avg) < eps
            out.ratio(ir) = NaN;
        else
            out.ratio(ir) = lift_avg / work_avg;
        end
        out.missing_steps(ir) = missing_steps;
    end
end
