%% Analyze Tropical Histogram Diagnostics (+4K)
% Thin wrapper: orchestration is in analysis/run_histogram_precip_work_analysis.m

clear all; close all; clc;

this_file = mfilename('fullpath');
if isempty(this_file)
    script_dir = pwd;
else
    script_dir = fileparts(this_file);
end

if isfolder(fullfile(script_dir, 'lib')) && isfolder(fullfile(script_dir, 'analysis'))
    matlab_root = script_dir;
else
    matlab_root = fileparts(script_dir);
end

addpath(fullfile(matlab_root, 'lib'));
addpath(fullfile(matlab_root, 'analysis'));

cfg = struct();
cfg.simulation_name = 'plus4k';
cfg.ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_PLUS_4K_CO2_1270ppmv/hist_2020010300_2022011800.nc';
cfg.time_weight_fn = @compute_time_weights_plus4k;
cfg.xmax = 2e-2;
cfg.title_prefix = '+4K';

run_histogram_precip_work_analysis(cfg);
