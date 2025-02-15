%% SUMMARY
% This script is 

%% PREPARE SOME NECESSARY CONSTANTS
clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_trial    = "trial_0006_Session2_02";

% declare some of the important paths
path_function   = fullfile(path_root, "functions");
path_outputs    = fullfile(path_root, "outputs");
path_bonestl    = fullfile(path_root, "data", "ct", "bone");
path_markerstl  = fullfile(path_root, "data", "ct", "marker");
addpath(genpath(path_function));

%% PREPARING FIGS

% prepare the figure object for visualizing 3d signal
fig1 = figure("Name", "3D Signal");
ax1  = axes("Parent", fig1);
grid(ax1, "on");
axis(ax1, "equal");
xlabel(ax1, "X");
ylabel(ax1, "Y");
zlabel(ax1, "Z");
view(ax1, [30 30]);
hold(ax1, "on");

%% READ MARKER

% get all stl files
file_markerSTLs     = dir(fullfile(path_markerstl, "*.stl"));
fullpath_markerSTLs = fullfile(path_markerstl, {file_markerSTLs.name});

% Initialize matrix to store centers
marker_centers = zeros(length(fullpath_markerSTLs), 3);

% Define color map for groups
colors = lines(4);
group_colors = containers.Map({
    'P_F_DIS', 'P_F_PRO', 'P_T_DIS', 'P_T_PRO'
}, num2cell(colors, 2));

markerpin = struct('name', "", 'centers', []);

% Process each marker file
for marker_idx = 1:length(fullpath_markerSTLs)

    % Determine group based on filename prefix
    [~, file_name, ~] = fileparts(fullpath_markerSTLs(marker_idx));
    prefix = regexpi(file_name, '^P_[FT]_(DIS|PRO)', 'match', 'once'); % Match valid prefixes

    % Validate the prefix
    if isempty(prefix) || ~isKey(group_colors, prefix)
        warning('Unexpected prefix "%s" found in file: %s', prefix, marker_files{marker_idx});
        continue;
    end


    % Read STL file
    markerSTL = stlread(fullpath_markerSTLs(marker_idx));
    % Extract points
    pointSTL = markerSTL.Points;
    % Fit a sphere to the points and calculate the center
    marker_centers(marker_idx, :) = estimateCenterMarker(pointSTL);
end


% Get the corresponding color
color = group_colors(prefix);
% scatter3(ax1, pointSTL(:, 1), pointSTL(:, 2), pointSTL(:, 3), 0.5, 'filled', 'MarkerFaceColor', color);
scatter3(ax1, marker_centers(:, 1), marker_centers(:, 2), marker_centers(:, 3), 20, 'filled', 'MarkerFaceColor', color);


