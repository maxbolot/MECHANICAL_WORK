function [time_indices, time_label] = select_time_indices(time_dt, temporal_mode, inst_date, range_start, range_end)
    % Selects time indices for instantaneous snapshots or inclusive date ranges.
    switch lower(string(temporal_mode))
        case "instantaneous"
            [~, idx] = min(abs(time_dt - inst_date));
            time_indices = idx;
            time_label = sprintf('instantaneous %s', datestr(time_dt(idx), 'yyyy-mm-dd HH:MM:SS UTC'));
        case "timerange"
            time_indices = find(time_dt >= range_start & time_dt <= range_end);
            if isempty(time_indices)
                error('No timestamps found inside range [%s, %s].', string(range_start), string(range_end));
            end
            time_label = sprintf('time-mean %s to %s', datestr(time_dt(time_indices(1)), 'yyyy-mm-dd'), datestr(time_dt(time_indices(end)), 'yyyy-mm-dd'));
        otherwise
            error('Unsupported temporal_mode: %s', temporal_mode);
    end
end
