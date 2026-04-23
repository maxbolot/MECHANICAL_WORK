if ~exist('matlab_base_dir', 'var')
    this_file = mfilename('fullpath');
    script_dir = fileparts(this_file);
    matlab_base_dir = fileparts(script_dir);
end
addpath(fullfile(matlab_base_dir, 'lib'));
addpath(fullfile(matlab_base_dir, 'analysis'));
addpath(fullfile(matlab_base_dir, 'presets'));

cfg.run_mode = 'non_thresholded';
% Supported: non_thresholded | prate_thresholded

cfg.threshold_percentile = 0.5;
cfg.g = 9.81;

run_hp_precip_delta_lonlat_map_analysis(cfg);
