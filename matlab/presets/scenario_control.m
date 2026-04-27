function scenario = scenario_control()
    scenario.name = 'control';
    scenario.label = 'Control';
    scenario.standard_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180/work_2020010300_2022011200.nc';
    scenario.standard_precip_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180/precip_2020010300_2022011200.nc';
    scenario.thresholded_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_prate_thresholded/work_2020010300_2022011200.nc';
    scenario.thresholded_precip_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded/precip_2020010300_2022011200.nc';
    scenario.thresholded_by_lat_band_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_prate_thresholded_by_lat_band/work_2020010300_2022011200.nc';
    scenario.thresholded_by_lat_band_precip_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded_by_lat_band/precip_2020010300_2022011200.nc';
    scenario.time_weight_fn = @compute_time_weights_control;
end