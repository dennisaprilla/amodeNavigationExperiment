%% PREPARE SOME NECESSARY CONSTANTS
clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
path_ukf     = "D:\Documents\MATLAB\unscentedkalmanfilter_registration\functions\ukf";
dir_function = "functions";
dir_outputs  = "outputs";

% declare some of the important paths
path_function = fullfile(path_root, dir_function);
path_outputs  = fullfile(path_root, dir_outputs);
path_bonestl  = fullfile(path_root, "data", "ct", "bone");
addpath(genpath(path_function));
addpath(path_ukf);

% get the screen size
scr_size  = get(0, 'ScreenSize');
% prepare the figure object for visualizing 3d signal
fig1 = figure("Name", "3D Signal", "Position", [0 0 0.5*scr_size(3), scr_size(4)]);
ax1  = axes("Parent", fig1);
grid(ax1, "on");
axis(ax1, "equal");
xlabel(ax1, "X");
ylabel(ax1, "Y");
zlabel(ax1, "Z");
view(ax1, [30 30]);
hold(ax1, "on");

% color maps
colormap1 = parula(256);
colormap2 = hot(256);

% load the measurements
load(fullfile(path_outputs, "measurementoutput_withnav_2025-01-28_19-18-21.mat"));
n_group    = length(all_measurements);

allUS   = struct('peak_3d', [], 'peak_intensity', [], 'peak_direction', []);
femurUS = struct('peak_3d', [], 'peak_intensity', [], 'peak_direction', []);
tibiaUS = struct('peak_3d', [], 'peak_intensity', [], 'peak_direction', []);

% obtain the 3d point cloud
for group_idx=1:n_group
    % get the probes measuremnets
    all_probes = all_measurements(group_idx).probes;
    n_probes   = length(all_probes);

    currentgroup_peak3d = [];
    currentgroup_peakintensity = [];
    currentgroup_peakdirection = [];
    for probe_idx=1:n_probes
        % skip if the peak is empty
        if(isempty(all_probes(probe_idx).peak_3d_inref))
            continue;
        end
            
        % get the biggest peak, the first column
        tmp = all_probes(probe_idx).peak_3d_inref(:, 1);
        currentgroup_peak3d = [currentgroup_peak3d, tmp];
        % get the intensity of the biggest peak
        tmp = all_probes(probe_idx).peak_2d(:,1);
        currentgroup_peakintensity = [currentgroup_peakintensity, tmp];
        % get the "direction" of the points
        tmp1 = all_probes(probe_idx).peak_3d_inref(:, 1);
        tmp2 = all_probes(probe_idx).signal_3d_inref(:, 1);
        tmp3 = (tmp2-tmp1)/norm((tmp2-tmp1));
        currentgroup_peakdirection = [currentgroup_peakdirection, tmp3(1:3)];
    end

    allUS.peak_3d = [allUS.peak_3d, currentgroup_peak3d];
    allUS.peak_intensity = [allUS.peak_intensity, currentgroup_peakintensity];
    allUS.peak_direction = [allUS.peak_direction, currentgroup_peakdirection];

    strings = split(all_measurements(group_idx).groupname, '_');
    if (strcmp(strings{2}, 'F'))
        femurUS.peak_3d = [femurUS.peak_3d, currentgroup_peak3d];
        femurUS.peak_intensity = [femurUS.peak_intensity, currentgroup_peakintensity];
        femurUS.peak_direction = [femurUS.peak_direction, currentgroup_peakdirection];
    elseif(strcmp(strings{2}, 'T'))
        tibiaUS.peak_3d = [tibiaUS.peak_3d, currentgroup_peak3d];
        tibiaUS.peak_intensity = [tibiaUS.peak_intensity, currentgroup_peakintensity];
        tibiaUS.peak_direction = [tibiaUS.peak_direction, currentgroup_peakdirection];
    end
end


%% FEMUR

%
% Plot the peak3d
scattersize_all  = femurUS.peak_intensity(2,:);
scattercolor_all = interp1(linspace(min(scattersize_all), max(scattersize_all), size(colormap1, 1)), colormap1, scattersize_all);
scatter3( ax1, ...
          femurUS.peak_3d(1,:), ...
          femurUS.peak_3d(2,:), ...
          femurUS.peak_3d(3,:), ...
          30, scattersize_all, 'filled');
% quiver3( ax1, ...
%           femurUS.peak_3d(1,:), ...
%           femurUS.peak_3d(2,:), ...
%           femurUS.peak_3d(3,:), ...
%           femurUS.peak_direction(1,:), ...
%           femurUS.peak_direction(2,:), ...
%           femurUS.peak_direction(3,:), ...
%           0.25, 'b');

% Extract all group names from the structure and find indices where the 
% middle part is 'F'
groupnames = {all_bonescans.groupname};
indices = find(cellfun(@(x) strcmp(x{2}, 'F'), cellfun(@(x) strsplit(x, '_'), groupnames, 'UniformOutput', false)));
for group_idx=indices
    % display the points
    scatter3(  ax1, ...
               all_bonescans(group_idx).voxels(1,:), ...
               all_bonescans(group_idx).voxels(2,:),  ...
               all_bonescans(group_idx).voxels(3,:), ...
               1, 'r', ...
              'filled', ...
              'Tag', "bone_surface");
end

% load stl files
femurCTstl_ref     = stlread(fullfile(path_bonestl, "Femur_adaptiveremeshed.stl"));
femurCTpc_ref      = pointCloud(femurCTstl_ref.Points);
femurCTnormals_ref = STLVertexNormals(femurCTstl_ref.ConnectivityList, femurCTstl_ref.Points)';
femurCTpoints_ref  = [femurCTpc_ref.Location'; ones(1, length(femurCTpc_ref.Location))];
% scatter3( ax1, ...
%           femurCTpoints_ref(1,:), ...
%           femurCTpoints_ref(2,:), ...
%           femurCTpoints_ref(3,:), ...
%           1, "red", "filled", ...
%           "Tag", "femur_original");
% quiver3( ax1, ...
%          femurCTpoints_ref(1,:), ...
%          femurCTpoints_ref(2,:), ...
%          femurCTpoints_ref(3,:), ...
%          femurCTnormals_ref(1,:), ...
%          femurCTnormals_ref(2,:), ...
%          femurCTnormals_ref(3,:), ...
%           "Tag", "femur_original");

% PCA alignment -----------------------------------------------------------

% get the principal component of femur stl
coeff = pca(femurCTstl_ref.Points);
mu    = mean(femurCTstl_ref.Points);
T_pcaFemurCT_ref = [coeff, mu'; 0 0 0 1];
% display_axis(ax1, T_pcaFemurCT_ref(1:3,4), T_pcaFemurCT_ref(1:3,1:3), 100, "Tag", "baseaxis_femurCT")

% get the principal componen of femur measurement
R_correction = eul2rotm([0 pi 0], "ZYX");
% t_correction = [-40 55 -15];
t_correction = [-10 50 -20];
coeff = pca(femurUS.peak_3d(1:3,:)') * R_correction;
mu    = mean(femurUS.peak_3d(1:3,:)') + t_correction; 
T_pcaFemurUS_ref = [coeff, mu'; 0 0 0 1];
% display_axis(ax1, T_pcaFemurUS_ref(1:3,4), T_pcaFemurUS_ref(1:3,1:3), 100, "Tag", "baseaxis_femurUS");

T_pcaFemurUS_pcaFemurCT  = T_pcaFemurUS_ref * inv(T_pcaFemurCT_ref);
femurCTpoints_pcaFemurCT = T_pcaFemurUS_pcaFemurCT * femurCTpoints_ref;
% scatter3( ax1, ...
%           femurCTpoints_pcaFemurCT(1,:), ...
%           femurCTpoints_pcaFemurCT(2,:), ...
%           femurCTpoints_pcaFemurCT(3,:), ...
%           0.5, "magenta", "filled", ...
%           "Tag", "femur_pca1");

% R-ICP alignment ---------------------------------------------------------
disp("Registering Femur...");

pc1 = pointCloud(femurUS.peak_3d(1:3,:)');
pc2 = pointCloud(femurCTpoints_pcaFemurCT(1:3, :)');

% pc1 and pc2 are defined somewhere above
iteration = 300;

% Preallocate for parfor
rmse_icp_all = zeros(1, iteration);
rmse_me_all  = zeros(1, iteration);
T_icp_all    = cell(1, iteration);
T_random_all = cell(1, iteration);

for i = 1:iteration
    try
        % make a random transformation
        T_random = rigidtform3d(randomRigidTransfrom(10,10));
        % apply icp
        [tform, movingReg, rmse_icp] = pcregistericp(pc2, pc1, ...
            "InitialTransform", T_random, ...
            "InlierRatio",      0.15, ...
            "MaxIterations",    50, ...
            "Verbose",          false);
        T_icp = tform.A;

        % search nearest neighbor 
        [nn_idx, nn_dist] = knnsearch(movingReg.Location, pc1.Location, "K", 20);
        rmse_me = mean(mean(nn_dist, 2),1);
        
        rmse_icp_all(i) = rmse_icp;
        rmse_me_all(i)  = rmse_me;
        T_icp_all{i}    = tform.A;
        T_random_all{i} = T_random.A;
    catch
        rmse_icp_all(i) = NaN;
        rmse_me_all(i)  = NaN;
        T_icp_all{i}    = [];
        T_random_all{i} = [];
    end
end

% search for the least rmse
[min_val, min_idx] = min(rmse_me_all);
T_icpFemur_pcaFemurUS   = T_icp_all{min_idx};

femurCTpoints_pcaFemurUS  = T_icpFemur_pcaFemurUS * T_pcaFemurUS_pcaFemurCT * femurCTpoints_ref;
femurCTnormals_pcaFemurUS = T_icpFemur_pcaFemurUS(1:3, 1:3) * T_pcaFemurUS_pcaFemurCT(1:3, 1:3) * femurCTnormals_ref;
scatter3( ax1, ...
          femurCTpoints_pcaFemurUS(1,:), ...
          femurCTpoints_pcaFemurUS(2,:), ...
          femurCTpoints_pcaFemurUS(3,:), ...
          0.5, "green", "filled", ...
          "Tag", "femur_ricp");
% quiver3( ax1, ...
%          femurCTpoints_pcaFemurUS(1,:), ...
%          femurCTpoints_pcaFemurUS(2,:), ...
%          femurCTpoints_pcaFemurUS(3,:), ...
%          femurCTnormals_pcaFemurUS(1,:), ...
%          femurCTnormals_pcaFemurUS(2,:), ...
%          femurCTnormals_pcaFemurUS(3,:), ...
%           "Tag", "femur_original");

% fine alignment
[T_all, mean_dist, history] = ukf_isotropic_registration_ex2( femurUS.peak_3d(1:3,:), femurCTpoints_pcaFemurUS(1:3,:), ...
                                   'movingnormal', femurUS.peak_direction, ...
                                   'fixednormal', femurCTnormals_pcaFemurUS, ...
                                   'normalratio', 0.5, ...
                                   'threshold', 0.0001, ...
                                   'iteration', 300, ...
                                   'expectednoise', 15, ...
                                   'sigmaxanneal', 0.99, ...
                                   'sigmaxtrans', 10, ...
                                   'sigmaxtheta', 10, ...
                                   'bestrmse', true, ...
                                   'verbose', false, ...
                                   'display', false);

% inverse the transformation (because, for some unknown reason, i made the
% registration from US measurement to the model)
T_ukfFemur_icpFemur = inv(T_all);
% transforrm
femurCTpoints_final  = T_ukfFemur_icpFemur * T_icpFemur_pcaFemurUS* T_pcaFemurUS_pcaFemurCT * femurCTpoints_ref;
femurCTnormals_final = T_ukfFemur_icpFemur(1:3, 1:3) * T_icpFemur_pcaFemurUS(1:3, 1:3) * T_pcaFemurUS_pcaFemurCT(1:3, 1:3) * femurCTnormals_ref;
scatter3( ax1, ...
          femurCTpoints_final(1,:), ...
          femurCTpoints_final(2,:), ...
          femurCTpoints_final(3,:), ...
          0.5, "blue", "filled", ...
          "Tag", "femur_ricp");
%


%% TIBIA

%
% plot the peak3d
scattersize_all  = tibiaUS.peak_intensity(2,:);
scattercolor_all = interp1(linspace(min(scattersize_all), max(scattersize_all), size(colormap1, 1)), colormap1, scattersize_all);
scatter3( ax1, ...
          tibiaUS.peak_3d(1,:), ...
          tibiaUS.peak_3d(2,:), ...
          tibiaUS.peak_3d(3,:), ...
          30, scattersize_all, 'filled');
colormap(ax1, "hot");
colorbar(ax1);

% Extract all group names from the structure and find indices where the 
% middle part is 'T'
groupnames = {all_bonescans.groupname};
indices = find(cellfun(@(x) strcmp(x{2}, 'T'), cellfun(@(x) strsplit(x, '_'), groupnames, 'UniformOutput', false)));
for group_idx=indices
    % display the points
    scatter3(  ax1, ...
               all_bonescans(group_idx).voxels(1,:), ...
               all_bonescans(group_idx).voxels(2,:),  ...
               all_bonescans(group_idx).voxels(3,:), ...
               1, 'r', ...
              'filled', ...
              'Tag', "bone_surface");
end

% load stl files
tibiaCTstl_ref     = stlread(fullfile(path_bonestl, "Tibia_adaptiveremeshed.stl"));
tibiaCTpc_ref      = pointCloud(tibiaCTstl_ref.Points);
tibiaCTnormals_ref = STLVertexNormals(tibiaCTstl_ref.ConnectivityList, tibiaCTstl_ref.Points)';
tibiaCTpoints_ref  = [tibiaCTpc_ref.Location'; ones(1, length(tibiaCTpc_ref.Location))];
% scatter3( ax1, ...
%           tibiaCTpoints_ref(1,:), ...
%           tibiaCTpoints_ref(2,:), ...
%           tibiaCTpoints_ref(3,:), ...
%           1, "red", "filled");

% PCA alignment -----------------------------------------------------------

% get the principal component of tibia stl
coeff = pca(tibiaCTstl_ref.Points);
mu    = mean(tibiaCTstl_ref.Points);
T_pcaTibiaCT_ref = [coeff, mu'; 0 0 0 1];   
% display_axis(ax1, T_pcaTibiaCT_ref(1:3,4), T_pcaTibiaCT_ref(1:3,1:3), 100, "Tag", "baseaxis_tibiaCT")

% get the principal componen of tibia measurement
R_correction = eul2rotm([0 pi 0], "ZYX");
% t_correction = [30 -20 -30];
t_correction = [-20 -10 -40];
coeff = pca(tibiaUS.peak_3d(1:3,:)') * R_correction;
mu    = mean(tibiaUS.peak_3d(1:3,:)') + t_correction; 
T_pcaTibiaUS_ref = [coeff, mu'; 0 0 0 1];
% display_axis(ax1, T_pcaTibiaUS_ref(1:3,4), T_pcaTibiaUS_ref(1:3,1:3), 100, "Tag", "baseaxis_tibiaUS")

T_pcaTibiaUS_pcaTibiaCT  = T_pcaTibiaUS_ref * inv(T_pcaTibiaCT_ref);
tibiaCTpoints_pcaTibiaCT = T_pcaTibiaUS_pcaTibiaCT * tibiaCTpoints_ref;
% scatter3( ax1, ...
%           tibiaCTpoints_pcaTibiaCT(1,:), ...
%           tibiaCTpoints_pcaTibiaCT(2,:), ...
%           tibiaCTpoints_pcaTibiaCT(3,:), ...
%           1, "magenta", "filled");


% R-ICP alignment ---------------------------------------------------------
disp("Registering Tibia...");

pc1 = pointCloud(tibiaUS.peak_3d(1:3,:)');
pc2 = pointCloud(tibiaCTpoints_pcaTibiaCT(1:3, :)');

% pc1 and pc2 are defined somewhere above
iteration = 300;

% Preallocate for parfor
rmse_icp_all = zeros(1, iteration);
rmse_me_all  = zeros(1, iteration);
T_icp_all    = cell(1, iteration);
T_random_all = cell(1, iteration);

for i = 1:iteration
    try
        % make a random transformation
        T_random = rigidtform3d(randomRigidTransfrom(10,10));
        % apply icp
        [tform, movingReg, rmse_icp] = pcregistericp(pc2, pc1, ...
            "InitialTransform", T_random, ...
            "InlierRatio",      0.20, ...
            "MaxIterations",    50, ...
            "Verbose",          false);
        T_icp = tform.A;

        % search nearest neighbor 
        [nn_idx, nn_dist] = knnsearch(movingReg.Location, pc1.Location, "K", 20);
        rmse_me = mean(mean(nn_dist, 2),1);
        
        rmse_icp_all(i) = rmse_icp;
        rmse_me_all(i)  = rmse_me;
        T_icp_all{i}    = tform.A;
        T_random_all{i} = T_random.A;
    catch
        rmse_icp_all(i) = NaN;
        rmse_me_all(i)  = NaN;
        T_icp_all{i}    = [];
        T_random_all{i} = [];
    end
end

% search for the least rmse
[min_val, min_idx] = min(rmse_me_all);
T_icpTibia_pcaTibiaUS   = T_icp_all{min_idx};

tibiaCTpoints_pcaTibiaUS = T_icpTibia_pcaTibiaUS * T_pcaTibiaUS_pcaTibiaCT * tibiaCTpoints_ref;
tibiaCTnormals_pcaTibiaUS = T_icpTibia_pcaTibiaUS(1:3, 1:3) * T_pcaTibiaUS_pcaTibiaCT(1:3, 1:3) * tibiaCTnormals_ref;
scatter3( ax1, ...
          tibiaCTpoints_pcaTibiaUS(1,:), ...
          tibiaCTpoints_pcaTibiaUS(2,:), ...
          tibiaCTpoints_pcaTibiaUS(3,:), ...
          1, "green", "filled", ...
          "Tag", "femur_ricp");
% quiver3( ax1, ...
%          tibiaCTpoints_pcaTibiaUS(1,:), ...
%          tibiaCTpoints_pcaTibiaUS(2,:), ...
%          tibiaCTpoints_pcaTibiaUS(3,:), ...
%          tibiaCTnormals_pcaTibiaUS(1,:), ...
%          tibiaCTnormals_pcaTibiaUS(2,:), ...
%          tibiaCTnormals_pcaTibiaUS(3,:), ...
%           "Tag", "femur_original")

% fine alignment
% show
[T_all, mean_dist, history] = ukf_isotropic_registration_ex2( tibiaUS.peak_3d(1:3,:), tibiaCTpoints_pcaTibiaUS(1:3,:), ...
                                   'movingnormal', tibiaUS.peak_direction, ...
                                   'fixednormal', tibiaCTnormals_pcaTibiaUS, ...
                                   'normalratio', 0.5, ...
                                   'threshold', 0.0001, ...
                                   'iteration', 300, ...
                                   'expectednoise', 15, ...
                                   'sigmaxanneal', 0.99, ...
                                   'sigmaxtrans', 10, ...
                                   'sigmaxtheta', 10, ...
                                   'bestrmse', true, ...
                                   'verbose', false, ...
                                   'display', false);

% inverse the transformation (because, for some unknown reason, i made the
% registration from US measurement to the model)
T_ukfTibia_icpTibia = inv(T_all);
% transforrm
tibiaCTpoints_final  = T_ukfTibia_icpTibia * T_icpTibia_pcaTibiaUS * T_pcaTibiaUS_pcaTibiaCT * tibiaCTpoints_ref;
tibiaCTnormals_final = T_ukfTibia_icpTibia(1:3, 1:3) * T_icpTibia_pcaTibiaUS(1:3, 1:3) * T_pcaTibiaUS_pcaTibiaCT(1:3, 1:3) * tibiaCTnormals_ref;
scatter3( ax1, ...
          tibiaCTpoints_final(1,:), ...
          tibiaCTpoints_final(2,:), ...
          tibiaCTpoints_final(3,:), ...
          1, "blue", "filled", ...
          "Tag", "femur_ricp");
%




