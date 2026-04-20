if ~exist('matlab_base_dir', 'var')
    this_file = mfilename('fullpath');
    script_dir = fileparts(this_file);
    matlab_base_dir = fileparts(script_dir);
end
addpath(fullfile(matlab_base_dir, 'lib'));
addpath(fullfile(matlab_base_dir, 'analysis'));
addpath(fullfile(matlab_base_dir, 'presets'));

cfg.scenario = scenario_control();
cfg.region = region_global();
cfg.smooth_before_plotting = true;
cfg.smooth_window_days = 5;
run_thresholded_region_analysis(cfg);