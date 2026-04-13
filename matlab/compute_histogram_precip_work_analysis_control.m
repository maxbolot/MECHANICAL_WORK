%% Analyze Tropical Histogram Diagnostics (Control)
% Reads tropical histogram diagnostics and computes:
%   1) Time-averaged area-weighted precipitation frequency.
%   2) Time-averaged precipitation-binned lift/work contributions from
%      normalized joint (precip, work) histograms.

clear all; close all; clc;

% Add local helper library folder to MATLAB path.
this_file = mfilename('fullpath');
if isempty(this_file)
    script_dir = pwd;
else
    script_dir = fileparts(this_file);
end
addpath(fullfile(script_dir, 'lib'));

% Input histogram file (tropical domain histogram output)
ncfile = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms/hist_2020010300_2022011200.nc';

if ~isfile(ncfile)
    error('File not found: %s', ncfile);
end

fprintf('Reading NetCDF file: %s\n', ncfile);

% Read coordinates and histogram fields
try
    time = ncread(ncfile, 'time');
    pr_edges = double(ncread(ncfile, 'pr_edges'));
    work_edges = double(ncread(ncfile, 'work_edges'));

    % ncread reverses NetCDF dimension order (Fortran/column-major):
    %   hist_area:  NetCDF (time, nbin_pr) -> MATLAB (nbin_pr, time)
    %   hist2d_*:   NetCDF (time, nbin_work, nbin_pr) -> MATLAB (nbin_pr, nbin_work, time)
    hist_area   = double(ncread(ncfile, 'hist_area'));
    hist2d_lift = double(ncread(ncfile, 'hist2d_lift'));
    hist2d_work = double(ncread(ncfile, 'hist2d_work'));

    fprintf('time size: %d\n', numel(time));
    fprintf('hist_area size (raw): %s\n', mat2str(size(hist_area)));
    fprintf('hist2d_lift size (raw): %s\n', mat2str(size(hist2d_lift)));
    fprintf('hist2d_work size (raw): %s\n', mat2str(size(hist2d_work)));

    % Permute to (time, ...) convention used throughout the script.
    hist_area   = hist_area';                          % (time, nbin_pr)
    hist2d_lift = permute(hist2d_lift, [3, 2, 1]);     % (time, nbin_work, nbin_pr)
    hist2d_work = permute(hist2d_work, [3, 2, 1]);     % (time, nbin_work, nbin_pr)

    fprintf('hist_area size (permuted): %s\n', mat2str(size(hist_area)));
    fprintf('hist2d_lift size (permuted): %s\n', mat2str(size(hist2d_lift)));
catch ME
    error('Failed to read required histogram variables: %s', ME.message);
end

ntime = numel(time);

% Build bin centers from edges.
% pr edges are linear in log-coordinate, so geometric centers are natural.
pr_centers = sqrt(pr_edges(1:end-1) .* pr_edges(2:end));
work_centers = 0.5 * (work_edges(1:end-1) + work_edges(2:end));

%% 1) Area-weighted precipitation frequency
% Normalize each time step so sum over precipitation bins is 1.
hist_area_sum = sum(hist_area, 2, 'omitnan');
hist_area_freq = hist_area ./ hist_area_sum;
hist_area_freq(hist_area_sum <= 0, :) = NaN;

% Schedule-aware time weighting (same logic as compute_global_work_average_control.m)
[time_weights_days, missing_steps] = compute_time_weights_control(time, ncfile);
if numel(time_weights_days) ~= ntime
    error('Time-weight vector length (%d) does not match data time length (%d).', numel(time_weights_days), ntime);
end

pr_freq_timeavg = weighted_time_mean_rows(hist_area_freq, time_weights_days);

%% 2) 1st moment of joint histograms -> precipitation-binned mean lift/work
% The conditional mean in precipitation bin pr is:
%   mean_lift(pr) = sum_w [ hist2d_lift(t,w,pr) * work_center(w) ] / sum_w [ hist2d_lift(t,w,pr) ]
%
% Forming this ratio per time step and then averaging is noisy at high-pr bins
% where only a handful of grid cells contribute and the per-step denominator is
% near-zero. Instead, accumulate numerator and denominator separately with time
% weights, then divide once (population-mean estimator).
work_centers_bcast = reshape(work_centers, 1, [], 1);  % (1, nbin_work, 1)

lift_numerator   = squeeze(sum(hist2d_lift .* work_centers_bcast, 2, 'omitnan'));  % (time, nbin_pr)
lift_denominator = squeeze(sum(hist2d_lift,                        2, 'omitnan'));  % (time, nbin_pr)
work_numerator   = squeeze(sum(hist2d_work .* work_centers_bcast, 2, 'omitnan'));  % (time, nbin_pr)
work_denominator = squeeze(sum(hist2d_work,                        2, 'omitnan'));  % (time, nbin_pr)

% Weighted accumulation: zero out timesteps where a bin has no area so they
% contribute nothing to either accumulator.
w_col = time_weights_days(:);  % (time, 1)

lift_num_w = lift_numerator   .* w_col;  lift_num_w(lift_denominator <= 0) = 0;
lift_den_w = lift_denominator .* w_col;  lift_den_w(lift_denominator <= 0) = 0;
work_num_w = work_numerator   .* w_col;  work_num_w(work_denominator <= 0) = 0;
work_den_w = work_denominator .* w_col;  work_den_w(work_denominator <= 0) = 0;

lift_by_pr_timeavg = sum(lift_num_w, 1) ./ sum(lift_den_w, 1);
work_by_pr_timeavg = sum(work_num_w, 1) ./ sum(work_den_w, 1);
lift_by_pr_timeavg(sum(lift_den_w, 1) <= 0) = NaN;
work_by_pr_timeavg(sum(work_den_w, 1) <= 0) = NaN;

%% Display summary
fprintf('\n=== RESULTS ===\n');
fprintf('Simulation: control\n');
fprintf('Time weighting: schedule-aware control (5-day before 2021-05-27, 1-day on/after)\n');
fprintf('Detected missing timesteps: %d\n', missing_steps);
fprintf('Number of time steps: %d\n', ntime);

%% Plot outputs
figure('Position', [120, 120, 1000, 700]);
set(gcf, 'color', 'w');
set(gcf, 'WindowStyle', 'docked');

subplot(2,1,1);
pr_freq_plot = pr_freq_timeavg;
pr_freq_plot(pr_freq_plot <= 0) = NaN;
plot(pr_centers, pr_freq_plot, 'k-', 'LineWidth', 1.8);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([min(pr_centers), 1.5e-2]);
xlabel('Precipitation bin center');
ylabel('Frequency');
title('Control: Time-Mean Area-Weighted Precipitation Frequency');
grid on;

subplot(2,1,2);
lift_plot = lift_by_pr_timeavg;
work_plot = work_by_pr_timeavg;
lift_plot(lift_plot <= 0) = NaN;
work_plot(work_plot <= 0) = NaN;
ratio_plot = lift_by_pr_timeavg ./ work_by_pr_timeavg;
ratio_plot(~isfinite(ratio_plot) | ratio_plot <= 0) = NaN;

yyaxis left;
h1 = plot(pr_centers, lift_plot, 'r-', 'LineWidth', 1.8); hold on;
h2 = plot(pr_centers, work_plot, 'b-', 'LineWidth', 1.8);
set(gca, 'XScale', 'log', 'YScale', 'log');
ylabel('Power density (W/m^2)');

yyaxis right;
h3 = plot(pr_centers, ratio_plot, 'k--', 'LineWidth', 1.5);
ylabel('Lift/Work ratio');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

xlim([min(pr_centers), 1.5e-2]);
xlabel('Precipitation bin center');
title('Control: Time-Mean Lift and Work in Precipitation Bins');
legend([h1 h2 h3], {'Lift', 'Work', 'Lift/Work'}, 'Location', 'best');
grid on;
hold off;

function row_timeavg = weighted_time_mean_rows(data_time_by_bin, time_weights)
    % Compute weighted mean over time independently for each bin/column.
    % Inputs:
    %   data_time_by_bin: (time, nbins)
    %   time_weights:     (time, 1)

    data_time_by_bin = double(data_time_by_bin);
    time_weights = double(time_weights(:));

    if size(data_time_by_bin, 1) ~= numel(time_weights)
        error('Size mismatch: data has %d timesteps but weights has %d.', size(data_time_by_bin, 1), numel(time_weights));
    end

    valid = ~isnan(data_time_by_bin) & isfinite(data_time_by_bin) & (time_weights > 0);
    weighted_data = data_time_by_bin .* time_weights;
    weighted_data(~valid) = 0;

    sum_w = sum(time_weights .* valid, 1);
    sum_xw = sum(weighted_data, 1);

    row_timeavg = sum_xw ./ sum_w;
    row_timeavg(sum_w <= 0) = NaN;
end
