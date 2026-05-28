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

cfg_ctrl.scenario = scenario_control();
cfg_ctrl.region   = region_tropics();

cfg_warm.scenario = scenario_plus4k();
cfg_warm.region   = region_tropics();

fprintf('=== Computing bin stats: Control ===\n');
stats_ctrl = compute_thresholded_bin_stats(cfg_ctrl);

fprintf('=== Computing bin stats: +4K ===\n');
stats_warm = compute_thresholded_bin_stats(cfg_warm);

plot_percentile_bin_ratio_comparison(stats_ctrl, stats_warm, cfg_ctrl.region.title_label, ...
    'Control', '+4K');
