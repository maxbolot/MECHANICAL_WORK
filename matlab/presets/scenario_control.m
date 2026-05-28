function scenario = scenario_control()
    data_root = get_globalfv3_data_root();

    scenario.name = 'control';
    scenario.label = 'Control';
    scenario.standard_ncfile = fullfile(data_root, 'work_coarse_C3072_360x180', 'work_2020010300_2022011200.nc');
    scenario.standard_precip_ncfile = fullfile(data_root, 'precip_coarse_C3072_360x180', 'precip_2020010300_2022011200.nc');
    scenario.thresholded_ncfile = fullfile(data_root, 'work_coarse_C3072_360x180_prate_thresholded', 'work_2020010300_2022011200.nc');
    scenario.thresholded_precip_ncfile = fullfile(data_root, 'precip_coarse_C3072_360x180_prate_thresholded', 'precip_2020010300_2022011200.nc');
    scenario.thresholded_by_lat_band_ncfile = fullfile(data_root, 'work_coarse_C3072_360x180_prate_thresholded_by_lat_band', 'work_2020010300_2022011200.nc');
    scenario.thresholded_by_lat_band_precip_ncfile = fullfile(data_root, 'precip_coarse_C3072_360x180_prate_thresholded_by_lat_band', 'precip_2020010300_2022011200.nc');
    scenario.time_weight_fn = @compute_time_weights_control;
end