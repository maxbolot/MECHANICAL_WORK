function values = read_required_var(ncfile, varname)
    % Reads a required variable and reports available alternatives on failure.
    try
        values = ncread(ncfile, varname);
    catch ME
        error('Error reading %s variable: %s\nAvailable variables: %s', varname, ME.message, get_available_vars(ncfile));
    end
end
