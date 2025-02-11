clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_function = "functions";
dir_trial    = "trial_0006_Session2_02";

% declare some of the important paths
path_function = fullfile(path_root, dir_function);
path_snapshot = fullfile(path_root, "data", dir_trial, "snapshot");
path_bonescan = fullfile(path_root, "data", dir_trial, "bonescan");
addpath(genpath(path_function));

%%


filename_bonescan = "VolumeOutput_2024-12-18_14-00-35.mha";
reader = MHAReader(fullfile(path_bonescan, filename_bonescan));
if reader.readVolumeImage()
    headerInfo = reader.getMHAHeader();
    volumeData = reader.getMHAVolume();
end

% find indices over threshold
threshold = 100;
idx_overthresh = find(volumeData>threshold);

% convert indices into voxel
voxelcoordinate = ind2sub_c_style(idx_overthresh, headerInfo.DimSize);
voxelcoordinate_homogeneous = [voxelcoordinate'; ones(1, length(voxelcoordinate))];

% get the transformation value from MHA reader
init_t = headerInfo.Offset;
init_s = diag(headerInfo.ElementSpacing);
init_R = reshape(headerInfo.TransformMatrix, [3,3]);
init_A = [init_R.*init_s, init_t'; 0 0 0 1];

% trasform the points
voxelcoordinate_homogeneous = init_A * voxelcoordinate_homogeneous;

% display the points
fig1 = figure("Name", "Rigid Bodies");
ax1  = axes(fig1);
grid(ax1, "on");
axis(ax1, "equal");
title(ax1, "Bone Surface");
hold(ax1, "on");
scatter3(  ax1, ...
           voxelcoordinate_homogeneous(1,:), ...
           voxelcoordinate_homogeneous(2,:),  ...
           voxelcoordinate_homogeneous(3,:), ...
           10, ...
          'filled', ...
          'Tag', "bone_surface");

scatter3(  ax1, ...
           0, 0, 0, ...
           30, ...
          'filled', ...
          'Tag', "bone_surface");



