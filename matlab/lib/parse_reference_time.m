function ref_time = parse_reference_time(ref_str)
    formats = {'yyyy-MM-dd HH:mm:ss', 'yyyy-MM-dd HH:mm', 'yyyy-MM-dd', ...
               'yyyy-MM-dd''T''HH:mm:ss', 'yyyy-MM-dd''T''HH:mm'};

    for i = 1:numel(formats)
        try
            ref_time = datetime(ref_str, 'InputFormat', formats{i}, 'TimeZone', 'UTC');
            return;
        catch
        end
    end

    % Let MATLAB attempt automatic parsing as a last resort.
    ref_time = datetime(ref_str, 'TimeZone', 'UTC');
end
