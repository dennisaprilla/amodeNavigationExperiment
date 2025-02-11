clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_function = "functions";
dir_trial    = "trial_0003_Session1_02";

% declare some of the important paths
path_function = fullfile(path_root, dir_function);
path_snapshot = fullfile(path_root, "data", dir_trial, "snapshot");
addpath(genpath(path_function));

% declare some of important names
groups      = {"A_F_GTR", "A_F_MID", "A_F_LEP", "A_F_MEP", "A_T_LEP", "A_T_MEP", "A_T_MID", "A_T_MAL"};
color      = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};

% create base rotation to trasnform qualisys base vector to matlab
% Qualisys has y direction as up, MATLAB has z direction as up
R_tmp = eul2rotm([0 0 pi/2], "ZYX");
t_tmp = [0 0 0]';
baseRotation_Qualisys2Matlab = [R_tmp, t_tmp; 0 0 0 1];
    
% prepare the figure
fig1 = figure("Name", "Rigid Bodies");
ax1  = axes(fig1);
grid(ax1, "on");
axis(ax1, "equal");
title(ax1, "Rigid Bodies");
view(ax1, [45,15]);
hold(ax1, "on");

% for all groups...
for group_idx=1:length(groups)
    % get the current group
    current_group = groups{group_idx};
    % create the current directory
    directory = fullfile(path_snapshot, current_group);
    % read the rigid body data from the CSV files
    current_groupdata = readCSV_stylusSnapshot(directory);
    
    % for all the snapshot inside the data struct...
    for data_idx=1:length(current_groupdata)
        % get the trasnsformation of the stylus
        T = baseRotation_Qualisys2Matlab * current_groupdata(data_idx).Stylus.T;
        R = T(1:3, 1:3);
        t = T(1:3, 4);
    
        % plot each snapshot
        plot3(ax1, t(1), t(2), t(3), 'ok', 'MarkerFaceColor', color{group_idx}, 'MarkerSize', 10, 'Tag', "test");
        quiver3( ax1, t(1), t(2), t(3), ...
                 R(1,3)*10, R(2,3)*10, R(3,3)*10, ...
                 'LineWidth', 2, 'Color', 'blue', 'Tag', "test");
    end
end