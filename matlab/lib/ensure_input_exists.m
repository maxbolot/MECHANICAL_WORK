function ensure_input_exists(ncfile)
    % Fails early when a required NetCDF input file is missing.
    if ~isfile(ncfile)
        error('File not found: %s', ncfile);
    end
end
