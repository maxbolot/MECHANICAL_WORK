function scenario = scenario_control()
    scenario.name = 'control';
    scenario.label = 'Control';
    scenario.standard_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180/work_2020010300_2022011200.nc';
    scenario.thresholded_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_prate_thresholded/work_2020010300_2022011200.nc';
    scenario.time_weight_fn = @compute_time_weights_control;
end