if ~exist('matlab_base_dir', 'var')
    this_file = mfilename('fullpath');
    script_dir = fileparts(this_file);
    matlab_base_dir = fileparts(script_dir);
end
addpath(fullfile(matlab_base_dir, 'lib'));
addpath(fullfile(matlab_base_dir, 'analysis'));
addpath(fullfile(matlab_base_dir, 'presets'));

% User switches
cfg.simulation = 'control';
% Supported: control | warming | control_prate_thresholded | warming_prate_thresholded

cfg.temporal_mode = 'instantaneous';
% Supported: instantaneous | timerange

% Preset times for display
cfg.instantaneous_date = datetime(2021, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
cfg.range_start = datetime(2020, 8, 1, 0, 0, 0, 'TimeZone', 'UTC');
cfg.range_end = datetime(2020, 8, 30, 23, 59, 59, 'TimeZone', 'UTC');

% Used only for thresholded simulations.
cfg.threshold_percentile = 0.5;

run_work_ratio_lonlat_map_analysis(cfg);