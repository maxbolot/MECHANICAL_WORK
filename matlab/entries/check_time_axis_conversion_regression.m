if ~exist('matlab_base_dir', 'var')
    this_file = mfilename('fullpath');
    script_dir = fileparts(this_file);
    matlab_base_dir = fileparts(script_dir);
end
addpath(fullfile(matlab_base_dir, 'lib'));

fprintf('Running time-axis conversion regression checks...\n');

% Test 1: CF units "days since" path.
nc1 = fullfile(tempdir, sprintf('time_regression_days_%d.nc', feature('getpid')));
if isfile(nc1)
    delete(nc1);
end
nccreate(nc1, 'time', 'Dimensions', {'time', 3});
ncwrite(nc1, 'time', [0; 1; 2]);
ncwriteatt(nc1, 'time', 'units', 'days since 2020-01-01 00:00:00');

t1 = double(ncread(nc1, 'time'));
days1 = normalize_time_axis_to_datenum(t1, nc1);
assert(numel(days1) == 3, 'Test 1 failed: unexpected output length.');
assert(all(abs(diff(days1) - 1.0) < 1.0e-10), 'Test 1 failed: day-step mismatch for CF days.');

[plot_time1, label1] = build_time_axis(nc1, t1);
assert(isdatetime(plot_time1), 'Test 1 failed: build_time_axis should return datetime for CF units.');
assert(contains(label1, 'days since', 'IgnoreCase', true), 'Test 1 failed: label should include units.');

% Test 2: CF units "hours since" path.
nc2 = fullfile(tempdir, sprintf('time_regression_hours_%d.nc', feature('getpid')));
if isfile(nc2)
    delete(nc2);
end
nccreate(nc2, 'time', 'Dimensions', {'time', 3});
ncwrite(nc2, 'time', [0; 24; 48]);
ncwriteatt(nc2, 'time', 'units', 'hours since 2020-02-01 00:00:00');

t2 = double(ncread(nc2, 'time'));
days2 = normalize_time_axis_to_datenum(t2, nc2);
assert(all(abs(diff(days2) - 1.0) < 1.0e-10), 'Test 2 failed: day-step mismatch for CF hours.');

% Test 3: numeric YYYYMMDD.f fallback path (no file metadata).
t3 = [20200101.0; 20200102.5; 20200104.0];
days3 = normalize_time_axis_to_datenum(t3, '');
assert(abs(days3(2) - days3(1) - 1.5) < 1.0e-10, 'Test 3 failed: YYYYMMDD.f conversion mismatch (step 1).');
assert(abs(days3(3) - days3(2) - 1.5) < 1.0e-10, 'Test 3 failed: YYYYMMDD.f conversion mismatch (step 2).');

delete(nc1);
delete(nc2);

fprintf('All time-axis conversion regression checks passed.\n');