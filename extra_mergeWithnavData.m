%% SUMMARY
% This script is for merging data produced from step1_withnav_all.m for
% LABEL2 (cleaned) and LABEL3 (window updated). They are from two different
% sets of directory; the data are the same, but the files that contain
% window definitions are different. 
% 
% The script step1_withnav_all.m was already been written. It will read 
% through one set of directory, and produce the necessary data 
% (all_measurement) for next step. However, later in the stage, i 
% reconsider to collect the necessary data (all_measurement) but wihtout
% window updated (that is LABEL2). I am too lazy to change the script such
% that it can go back and forth between two directories. So what i did, i
% just run the script 2 times and have to different all_measurement.
%
% But that is ugly. Almost all data in the two all_measurement are exact 
% same, exept the detected 2d and 3d peak. So this script will be used 
% for merging those two all_measurement into one. I will add the data from 
% LABEL3 to peak_2d and peak_3d in the all_measurement from LABEL2

%% PREPARE SOME NECESSARY CONSTANTS

clear; clc; close all;

% change this
path_root  = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_trial  = "trial_0003_Session1_02";

% declare some of the important paths
path_function = fullfile(path_root, "functions");
path_outputs  = fullfile(path_root, "outputs");
path_snapshot = fullfile(path_root, "data", dir_trial, "snapshot");
path_bonescan = fullfile(path_root, "data", dir_trial, "bonescan");
addpath(genpath(path_function));

% specify LABEL2 and LABEL3 data
file_LABEL2 = "output-measurement_withnav-cleaned_2025-02-14-11-22-56.mat";
file_LABEL3 = "output-measurement_withnav-windowupdated_2025-02-14-11-23-58.mat";

% save?
is_saving = true;

%% MERGING

% get the LABEL's 3 all_measurement and store it into a temporary variable
load(fullfile(path_outputs, file_LABEL3));
tmp_measurements = all_measurements;

% get the LABEL2's all_measurement
load(fullfile(path_outputs, file_LABEL2));
n_groups = length(all_measurements);

% loop through all of the groups
for group_idx=1:n_groups
    % get the probes measuremnets
    all_probes = all_measurements(group_idx).probes;
    tmp_probes = tmp_measurements(group_idx).probes;
    n_probes   = length(all_probes);

    % loop through all of the probes
    for probe_idx=1:n_probes
        all_measurements(group_idx).probes(probe_idx).peak_2d       = [all_measurements(group_idx).probes(probe_idx).peak_2d,       tmp_measurements(group_idx).probes(probe_idx).peak_2d];
        all_measurements(group_idx).probes(probe_idx).window_2d     = [all_measurements(group_idx).probes(probe_idx).window_2d,     tmp_measurements(group_idx).probes(probe_idx).window_2d];
        all_measurements(group_idx).probes(probe_idx).peak_3d_inref = [all_measurements(group_idx).probes(probe_idx).peak_3d_inref, tmp_measurements(group_idx).probes(probe_idx).peak_3d_inref];
    end
end


%% SAVING

% save the data
if(is_saving)
    currentTime = datestr(now, 'yyyy-mm-dd-HH-MM-SS');

    str_LABEL2       = split(file_LABEL2, '_');
    str_LABEL2a      = split(str_LABEL2(2), '-');
    str_LABEL3       = split(file_LABEL3, '_');
    str_LABEL3a      = split(str_LABEL3(2), '-');

    file_output = str_LABEL2(1) + '_' + ...
                  str_LABEL2a(1)+'-'+str_LABEL2a(2)+'-'+str_LABEL3a(2)+'_' + ...
                  currentTime+'.mat';
    save(fullfile(path_outputs, file_output), 'all_measurements', 'all_bonescans');
end























