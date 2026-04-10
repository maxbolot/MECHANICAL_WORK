%% Compute work/lift for all regions and both simulations; plot legacy axis-1 style bars
% Regions: global, tropics, maritime continent, northern midlatitudes
% Simulations: control and +4K CO2 1270 ppmv

clear all; close all; clc;

% Add local helper library folder to MATLAB path.
this_file = mfilename('fullpath');
if isempty(this_file)
    script_dir = pwd;
else
    script_dir = fileparts(this_file);
end
addpath(fullfile(script_dir, 'lib'));

% Input files
control_file = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180/work_2020010800_2021011200.nc';
warming_file = '/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/work_2020010300_2021011600.nc';

if ~isfile(control_file)
    error('Control file not found: %s', control_file);
end
if ~isfile(warming_file)
    error('Warming file not found: %s', warming_file);
end

regions = {
    struct('name', 'global', 'latmin', -90, 'latmax', 90, 'lonmin', NaN, 'lonmax', NaN), ...
    struct('name', 'tropics', 'latmin', -30, 'latmax', 30, 'lonmin', NaN, 'lonmax', NaN), ...
    struct('name', 'maritime continent', 'latmin', -15, 'latmax', 15, 'lonmin', 90, 'lonmax', 150), ...
    struct('name', 'north. midlatitudes', 'latmin', 30, 'latmax', 60, 'lonmin', NaN, 'lonmax', NaN)
};

fprintf('Reading control simulation: %s\n', control_file);
control = process_file(control_file, regions, false);

fprintf('Reading warming simulation: %s\n', warming_file);
warming = process_file(warming_file, regions, true);

% Build stacked-bar input [lift, ke] in legacy row order:
% [control zone1; warming zone1; control zone2; warming zone2; ...]
Y = zeros(8, 2);
for i = 1:4
    Y(2*i - 1, :) = [control.lift(i), control.ke(i)];
    Y(2*i, :)     = [warming.lift(i), warming.ke(i)];
end

% Add spacer row after each region pair to mimic legacy axis 1 layout.
newY = zeros(12, 2);
newY([1 2 4 5 7 8 10 11], :) = Y;

%% Plot (legacy axis-1 style)
fig = figure('units', 'inch', 'position', [0, 0, 8, 3.8]);
set(gcf, 'color', 'w');

ax1 = axes(fig);
set(ax1, 'FontSize', 13);
b1 = bar(newY, 'stacked', 'FaceColor', 'flat');
set(gca, 'XLim', [0 12], 'XTick', [1.5 4.5 7.5 10.5], ...
    'XTickLabel', {'global', 'tropics', 'maritime continent', 'north. midlatitudes'});
ylabel('W m^{-2}');
title('Legacy Axis 1: work = lift + ke');

% Colors: control/warming hues with lighter second stack component.
control_dark = [0.10, 0.35, 0.80];
control_light = [0.65, 0.80, 0.98];
warming_dark = [0.85, 0.20, 0.20];
warming_light = [0.98, 0.72, 0.72];
white_bar = [1.0, 1.0, 1.0];

c1 = [control_dark; warming_dark; white_bar; control_dark; warming_dark; white_bar; ...
      control_dark; warming_dark; white_bar; control_dark; warming_dark; white_bar];
c2 = [control_light; warming_light; white_bar; control_light; warming_light; white_bar; ...
      control_light; warming_light; white_bar; control_light; warming_light; white_bar];

b1(1).CData = c1;  % lift
b1(2).CData = c2;  % ke = work - lift

legend({'lift', 'ke = work - lift'}, 'Location', 'northwest');
legend('boxoff');
grid on;

%% Print summary table
fprintf('\n=== Spatial + Time Averaged Results ===\n');
fprintf('%-22s %-12s %12s %12s %12s %12s %10s\n', 'Region', 'Simulation', 'work', 'lift', 'ke', 'lift/work', 'missing');
for i = 1:4
    fprintf('%-22s %-12s %12.6f %12.6f %12.6f %12.6f %10d\n', regions{i}.name, 'control', ...
        control.work(i), control.lift(i), control.ke(i), control.ratio(i), control.missing_steps(i));
    fprintf('%-22s %-12s %12.6f %12.6f %12.6f %12.6f %10d\n', regions{i}.name, 'warming', ...
        warming.work(i), warming.lift(i), warming.ke(i), warming.ratio(i), warming.missing_steps(i));
end


function out = process_file(ncfile, regions, use_schedule_weights)
    lat = double(ncread(ncfile, 'lat'));
    lon = double(ncread(ncfile, 'lon'));
    time = double(ncread(ncfile, 'time'));

    [work3d, dim_names] = read_var_as_lon_lat_time(ncfile, 'work');
    [lift3d, ~] = read_var_as_lon_lat_time(ncfile, 'lift');

    fprintf('  %s work dims interpreted as [lon, lat, time] from source dims [%s]\n', ...
        ncfile, strjoin(dim_names, ', '));

    if use_schedule_weights
        [time_weights_days, missing_steps] = compute_time_weights_plus4k(time, ncfile);
    else
        % Control simulation: uniform weighting over available timesteps.
        time_weights_days = ones(numel(time), 1);
        missing_steps = 0;
    end

    nreg = numel(regions);
    out.work = zeros(nreg, 1);
    out.lift = zeros(nreg, 1);
    out.ke = zeros(nreg, 1);
    out.ratio = zeros(nreg, 1);
    out.missing_steps = zeros(nreg, 1);

    lon_wrapped = mod(lon(:), 360);

    for ir = 1:nreg
        r = regions{ir};

        lat_idx = find(lat >= r.latmin & lat <= r.latmax);
        if isnan(r.lonmin)
            lon_idx = 1:numel(lon);
        else
            lon_idx = find(lon_wrapped >= r.lonmin & lon_wrapped <= r.lonmax);
        end

        if isempty(lat_idx) || isempty(lon_idx)
            error('No grid cells found for region: %s', r.name);
        end

        lat_weights = cosd(lat(lat_idx));
        w3 = reshape(lat_weights(:), 1, [], 1);

        wr = work3d(lon_idx, lat_idx, :);
        lr = lift3d(lon_idx, lat_idx, :);

        w_num = sum(wr .* w3, [1, 2], 'omitnan');
        w_den = sum(double(~isnan(wr)) .* w3, [1, 2], 'omitnan');
        work_ts = squeeze(w_num ./ w_den);

        l_num = sum(lr .* w3, [1, 2], 'omitnan');
        l_den = sum(double(~isnan(lr)) .* w3, [1, 2], 'omitnan');
        lift_ts = squeeze(l_num ./ l_den);

        work_avg = weighted_nanmean(work_ts, time_weights_days);
        lift_avg = weighted_nanmean(lift_ts, time_weights_days);

        out.work(ir) = work_avg;
        out.lift(ir) = lift_avg;
        out.ke(ir) = work_avg - lift_avg;
        if abs(work_avg) < eps
            out.ratio(ir) = NaN;
        else
            out.ratio(ir) = lift_avg / work_avg;
        end
        out.missing_steps(ir) = missing_steps;
    end
end
