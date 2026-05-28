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
rehash path;

cfg = struct();
cfg.ncfile = fullfile(get_globalfv3_data_root(), 'work_histograms', 'hist_2020010300_2022011200.nc');
cfg.time_weight_fn = @compute_time_weights_control;
cfg.title_prefix = 'Control';
cfg.make_plot = true;
cfg.print_summary = true;

% Optional axis limits: set to [] for auto.
cfg.xlim = [];
cfg.ylim_hist_area = [0 1e13];
cfg.ylim_area_fraction = [0 0.02];

run_tropics_precip_histogram_time_mean_analysis(cfg);
