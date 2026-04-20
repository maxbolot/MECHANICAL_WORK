function scenario = scenario_plus4k()
    scenario.name = 'plus4k';
    scenario.label = '+4K';
    scenario.standard_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/work_2020010300_2022011800.nc';
    scenario.thresholded_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded/work_2020010300_2022011800.nc';
    scenario.time_weight_fn = @compute_time_weights_plus4k;
end