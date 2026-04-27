if ~exist('matlab_base_dir', 'var')
    this_file = mfilename('fullpath');
    script_dir = fileparts(this_file);
    matlab_base_dir = fileparts(script_dir);
end
addpath(fullfile(matlab_base_dir, 'lib'));
addpath(fullfile(matlab_base_dir, 'analysis'));
addpath(fullfile(matlab_base_dir, 'presets'));

cfg.run_mode = 'wet_day_only';
% Supported: non_thresholded | wet_day_only | prate_thresholded | prate_thresholded_by_lat_band

% Percentile coordinate value used only when run_mode is thresholded.
% Example: 0.5 for p50, 0.9 for p90, 0.99 for p99 (if present in file).
cfg.threshold_percentile = NaN;
if any(strcmpi(cfg.run_mode, {'prate_thresholded', 'prate_thresholded_by_lat_band'}))
    cfg.threshold_percentile = 0.5;
end

% Wet-day filter threshold (kg m^-2 s^-1) used only when run_mode='wet_day_only'.
% Override this value as needed; default is 1e-5.
cfg.wet_day_threshold = 1.0e-5;

cfg.g = 9.81;

% Optional manual color limits per panel; leave empty for auto symmetric limits.
cfg.clim_hp = [-2500, 2500];
cfg.clim_precip = [-1e-4, 1e-4];

% Optional manual contour levels per panel; leave empty for auto levels.
cfg.contour_levels_hp = [3500 4000 4500];
cfg.contour_levels_precip = [15 20 100 150 200] .* 1e-6;

% Toggle contour labels on/off.
cfg.show_contour_labels = false;

run_hp_precip_delta_lonlat_map_analysis(cfg);
