function stats = build_percentile_bin_stats(percentile, work_tail_avg, lift_tail_avg, work_total_avg, lift_total_avg)
    % Converts tail-mean percentile diagnostics into non-overlapping bin statistics.
    p = double(percentile(:));
    nperc = numel(p);
    if nperc < 1
        error('No percentile thresholds found in the thresholded file.');
    end
    if any(diff(p) <= 0)
        error('Percentile coordinate must be strictly increasing.');
    end
    if p(1) <= 0 || p(end) >= 1
        error('Percentile coordinate values must be strictly inside (0, 1).');
    end

    lower = [0.0; p];
    upper = [p; 1.0];
    mass = upper - lower;
    nbin = numel(mass);

    work_bin_avg = nan(nbin, 1);
    lift_bin_avg = nan(nbin, 1);

    tail_mass_first = 1.0 - p(1);
    work_bin_avg(1) = (work_total_avg - tail_mass_first * work_tail_avg(1)) / mass(1);
    lift_bin_avg(1) = (lift_total_avg - tail_mass_first * lift_tail_avg(1)) / mass(1);

    for ib = 2:nperc
        p_lo = p(ib - 1);
        p_hi = p(ib);
        tail_mass_lo = 1.0 - p_lo;
        tail_mass_hi = 1.0 - p_hi;
        work_bin_avg(ib) = (tail_mass_lo * work_tail_avg(ib - 1) - tail_mass_hi * work_tail_avg(ib)) / mass(ib);
        lift_bin_avg(ib) = (tail_mass_lo * lift_tail_avg(ib - 1) - tail_mass_hi * lift_tail_avg(ib)) / mass(ib);
    end

    work_bin_avg(end) = work_tail_avg(end);
    lift_bin_avg(end) = lift_tail_avg(end);

    ratio_bin = lift_bin_avg ./ work_bin_avg;
    ratio_bin(abs(work_bin_avg) < eps) = NaN;

    labels = cell(nbin, 1);
    labels{1} = sprintf('<P%s', format_percentile_value_local(100.0 * p(1)));
    for ib = 2:nperc
        labels{ib} = sprintf('P%s-P%s', format_percentile_value_local(100.0 * p(ib - 1)), format_percentile_value_local(100.0 * p(ib)));
    end
    labels{end} = sprintf('>P%s', format_percentile_value_local(100.0 * p(end)));

    stats = struct();
    stats.bin_lower = lower;
    stats.bin_upper = upper;
    stats.bin_mass = mass;
    stats.bin_center = 0.5 * (lower + upper);
    stats.work_bin_avg = work_bin_avg;
    stats.lift_bin_avg = lift_bin_avg;
    stats.ratio_bin = ratio_bin;
    stats.labels = labels;
end


function text_out = format_percentile_value_local(value_in)
    text_out = regexprep(sprintf('%.6f', value_in), '0+$', '');
    text_out = regexprep(text_out, '\.$', '');
end
