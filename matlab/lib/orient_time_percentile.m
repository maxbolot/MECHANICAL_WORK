function values = orient_time_percentile(values_in, ntime, nperc)
    % Ensures averaged arrays are oriented as [time, percentile].
    values = values_in;
    if size(values, 1) ~= ntime && size(values, 2) == ntime
        values = values.';
    end
    if size(values, 1) ~= ntime || size(values, 2) ~= nperc
        error('Unexpected averaged field shape: got=%s expected=[%d %d]', mat2str(size(values)), ntime, nperc);
    end
end
