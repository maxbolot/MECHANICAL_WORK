function [ncfile, simulation_label, is_thresholded] = resolve_work_ratio_input_file(simulation)
    % Resolves scenario keywords to input files and thresholded mode flags.
    switch lower(string(simulation))
        case "control"
            sc = scenario_control();
            ncfile = sc.standard_ncfile;
            simulation_label = 'Control';
            is_thresholded = false;
        case "warming"
            sc = scenario_plus4k();
            ncfile = sc.standard_ncfile;
            simulation_label = 'Warming (+4K)';
            is_thresholded = false;
        case "control_prate_thresholded"
            sc = scenario_control();
            ncfile = sc.thresholded_ncfile;
            simulation_label = 'Control';
            is_thresholded = true;
        case "warming_prate_thresholded"
            sc = scenario_plus4k();
            ncfile = sc.thresholded_ncfile;
            simulation_label = 'Warming (+4K)';
            is_thresholded = true;
        otherwise
            error('Unsupported simulation: %s', simulation);
    end

    ensure_input_exists(ncfile);
end
