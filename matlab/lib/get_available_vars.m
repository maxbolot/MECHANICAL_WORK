function vars = get_available_vars(ncfile)
    try
        ncinfo_struct = ncinfo(ncfile);
        vars = {ncinfo_struct.Variables.Name};
        vars = sprintf('%s ', vars{:});
    catch
        vars = 'Unable to retrieve variable list';
    end
end
