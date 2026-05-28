this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = which(mfilename());
end
if isempty(this_file)
    error('Could not resolve path for %s.', mfilename());
end
script_dir = fileparts(this_file);
matlab_base_dir = fileparts(script_dir);
addpath(fullfile(matlab_base_dir, 'lib'));
addpath(fullfile(matlab_base_dir, 'analysis'));
addpath(fullfile(matlab_base_dir, 'presets'));
rehash path;

cfg.scenario = scenario_control();
cfg.region = region_tropics();

report = regression_check_thresholded_bin_stats(cfg); %#ok<NASGU>
