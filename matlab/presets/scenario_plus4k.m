function scenario = scenario_plus4k()
    scenario.name = 'plus4k';
    scenario.label = '+4K';
    scenario.standard_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/work_2020010300_2022011800.nc';
    scenario.standard_precip_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/precip_2020010300_2022012000.nc';
    scenario.thresholded_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded/work_2020010300_2022011800.nc';
    scenario.thresholded_precip_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded/precip_2020010300_2022012000.nc';
    scenario.thresholded_by_lat_band_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded_by_lat_band/work_2020010300_2022011800.nc';
    scenario.thresholded_by_lat_band_precip_ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded_by_lat_band/precip_2020010300_2022011800.nc';
    scenario.time_weight_fn = @compute_time_weights_plus4k;
end