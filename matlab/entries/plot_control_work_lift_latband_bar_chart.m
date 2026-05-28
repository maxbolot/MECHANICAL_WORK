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

cfg = struct();
cfg.scenario = scenario_control();
cfg.threshold_percentile = 0.99999;
cfg.make_plots = true;
cfg.print_summary = true;
cfg.window_style = 'docked';

run_control_work_lift_latband_bar_chart_analysis(cfg);
