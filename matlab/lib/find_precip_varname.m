function precip_varname = find_precip_varname(ncfile)
    % Detects the precipitation variable name from common naming conventions.
    info = ncinfo(ncfile);
    all_names = lower(string({info.Variables.Name}));

    candidates = ["pr", "prate", "prmax", "pratesfc", "pratesfc_coarse", ...
        "precip", "precipitation", "pr_sfc", "rainrate"];

    precip_varname = "";
    for i = 1:numel(candidates)
        idx = find(all_names == candidates(i), 1);
        if ~isempty(idx)
            precip_varname = string(info.Variables(idx).Name);
            return;
        end
    end

    contains_candidates = ["prmax", "pratesfc", "prate", "precip", "rainrate"];
    for i = 1:numel(contains_candidates)
        idx = find(contains(all_names, contains_candidates(i)), 1);
        if ~isempty(idx)
            precip_varname = string(info.Variables(idx).Name);
            return;
        end
    end

    error('No precipitation-like variable found in %s. Available variables: %s', ...
        ncfile, strjoin(cellstr(all_names), ', '));
end
