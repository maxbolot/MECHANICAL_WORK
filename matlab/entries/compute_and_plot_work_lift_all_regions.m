%% Compute work/lift for all regions and both simulations; plot legacy axis-1 style bars
% Thin wrapper: orchestration is in analysis/run_work_lift_all_regions_analysis.m

clear all; close all; clc;

this_file = mfilename('fullpath');
if isempty(this_file)
    script_dir = pwd;
else
    script_dir = fileparts(this_file);
end

if isfolder(fullfile(script_dir, 'lib')) && isfolder(fullfile(script_dir, 'analysis'))
    matlab_root = script_dir;
else
    matlab_root = fileparts(script_dir);
end

addpath(fullfile(matlab_root, 'lib'));
addpath(fullfile(matlab_root, 'analysis'));

cfg = struct();
run_work_lift_all_regions_analysis(cfg);
