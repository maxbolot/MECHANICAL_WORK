if ~exist('matlab_base_dir', 'var')
    this_file = mfilename('fullpath');
    script_dir = fileparts(this_file);
    matlab_base_dir = fileparts(script_dir);
end
addpath(fullfile(matlab_base_dir, 'lib'));
addpath(fullfile(matlab_base_dir, 'analysis'));
addpath(fullfile(matlab_base_dir, 'presets'));

cfg.run_mode = 'prate_thresholded';
% Supported: non_thresholded | prate_thresholded

cfg.threshold_percentile = 0.5;
cfg.g = 9.81;

% Optional manual color limits per panel; leave empty for auto limits.
cfg.clim_hp = [0, 2e4];
cfg.clim_precip = [];

% Optional manual contour levels per panel; leave empty for auto levels.
cfg.contour_levels_hp = [];
cfg.contour_levels_precip = [1e-10];

run_hp_precip_control_lonlat_map_analysis(cfg);
