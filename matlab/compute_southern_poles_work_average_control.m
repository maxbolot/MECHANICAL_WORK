%% Compute Mechanical Work and Lift Work Averages
% This script reads a NetCDF file and computes the regional average
% (90S to 60S) of mechanical work and lift work, averaged over time.

clear all; close all; clc;

% Add local helper library folder to MATLAB path.
this_file = mfilename('fullpath');
if isempty(this_file)
    script_dir = pwd;
else
    script_dir = fileparts(this_file);
end
addpath(fullfile(script_dir, 'lib'));

% Input file
ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180/work_2020010300_2022011200.nc';

% Check if file exists
if ~isfile(ncfile)
    error('File not found: %s', ncfile);
end

fprintf('Reading NetCDF file: %s\n', ncfile);

% Read latitude and check available variables
try
    lat = ncread(ncfile, 'lat');
    lon = ncread(ncfile, 'lon');
    time = ncread(ncfile, 'time');
    fprintf('File dimensions: lat=%d, lon=%d, time=%d\n', length(lat), length(lon), length(time));
catch ME
    error('Error reading dimensions: %s', ME.message);
end

% Select southern polar region indices:
% latitude: 90S to 60S, all longitudes.
lat_idx = find(lat >= -90 & lat <= -60);
lon_idx = 1:length(lon);

if isempty(lat_idx)
    error('No grid points found for southern polar bounds.');
end
fprintf('Southern poles selected: lat [%g, %g], all longitudes\n', ...
    min(lat(lat_idx)), max(lat(lat_idx)));

% Latitude weights for area-weighted averaging on a lat-lon grid
lat_selected = double(lat(lat_idx));
lat_weights = cosd(lat_selected(:));
if all(lat_weights == 0)
    error('Latitude weights are all zero in the selected southern polar region.');
end

% Read mechanical work variable
try
    work = ncread(ncfile, 'work');
    fprintf('Work variable shape: %s\n', mat2str(size(work)));
catch ME
    error('Error reading work variable: %s\nAvailable variables: %s', ME.message, get_available_vars(ncfile));
end

% Read lift variable
try
    lift = ncread(ncfile, 'lift');
    fprintf('Lift variable shape: %s\n', mat2str(size(lift)));
catch ME
    error('Error reading lift variable: %s', ME.message);
end

% Extract regional data
% Assuming dimensions are (lon, lat, time)
work_region = work(lon_idx, lat_idx, :);
lift_region = lift(lon_idx, lat_idx, :);

% Compute area-weighted spatial mean (cos(lat)) at each time step,
% then average in time.
weight_3d = reshape(lat_weights, 1, [], 1);

work_num = sum(work_region .* weight_3d, [1, 2], 'omitnan');
work_den = sum(double(~isnan(work_region)) .* weight_3d, [1, 2], 'omitnan');
work_spatial_avg = squeeze(work_num ./ work_den);

lift_num = sum(lift_region .* weight_3d, [1, 2], 'omitnan');
lift_den = sum(double(~isnan(lift_region)) .* weight_3d, [1, 2], 'omitnan');
lift_spatial_avg = squeeze(lift_num ./ lift_den);

[time_weights_days, missing_steps] = compute_time_weights_control(time, ncfile);
work_avg = weighted_nanmean(work_spatial_avg, time_weights_days);
lift_avg = weighted_nanmean(lift_spatial_avg, time_weights_days);

% Compute lift/work ratio for each time step and for the time-mean values
ratio_spatial_avg = lift_spatial_avg ./ work_spatial_avg;
ratio_spatial_avg(abs(work_spatial_avg) < eps) = NaN;

if abs(work_avg) < eps
    ratio_avg = NaN;
else
    ratio_avg = lift_avg / work_avg;
end

% Build plotting axis from NetCDF time variable (prefer datetime if units exist)
[plot_time, plot_time_label] = build_time_axis(ncfile, time);

% Display results
fprintf('\n');
fprintf('=== RESULTS ===\n');
fprintf('Region: Southern Poles (90°S to 60°S, all longitudes)\n');
fprintf('Spatial weighting: cos(latitude)\n');
fprintf('Time weighting: schedule-aware control (5-day before 2021-05-27, 1-day on/after)\n');
fprintf('Detected missing timesteps: %d\n', missing_steps);
fprintf('Number of time steps: %d\n', length(time));
fprintf('\nMechanical Work (averaged over southern poles and time): %.6f\n', work_avg);
fprintf('Lift Work (averaged over southern poles and time): %.6f\n', lift_avg);
fprintf('Lift/Work ratio (time-mean values): %.6f\n', ratio_avg);

% Create figure showing time series of regional means and ratio
figure('Position', [100, 100, 1000, 600]);
set(gcf, 'color', 'w');
set(gcf, 'WindowStyle', 'docked');

subplot(3, 1, 1);
plot(plot_time, work_spatial_avg, 'b-', 'LineWidth', 1.5);
xlabel(plot_time_label); ylabel('Work'); title('Mechanical Work - Southern Poles Average Time Series');
grid on;

subplot(3, 1, 2);
plot(plot_time, lift_spatial_avg, 'r-', 'LineWidth', 1.5);
xlabel(plot_time_label); ylabel('Lift'); title('Lift Work - Southern Poles Average Time Series');
grid on;

subplot(3, 1, 3);
plot(plot_time, ratio_spatial_avg, 'k-', 'LineWidth', 1.5);
xlabel(plot_time_label); ylabel('Lift/Work'); title('Lift/Work Ratio - Southern Poles Average Time Series');
ylim([0 1])
grid on;

% Save figure
% savefig('./southern_poles_work_analysis.fig');
% fprintf('\nFigure saved as: southern_poles_work_analysis.fig\n');
