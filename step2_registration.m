%% PREPARE SOME NECESSARY CONSTANTS
clear; clc; close all;

% [edit] some necessary path
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
path_ukf     = "D:\Documents\MATLAB\unscentedkalmanfilter_registration\functions\ukf";

% [edit] declare some of the important paths
path_function = fullfile(path_root, "functions");
path_outputs  = fullfile(path_root, "outputs");
path_bonestl  = fullfile(path_root, "data", "ct", "bone");

% add paths
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

%% OBTAINING 3D POINT CLOUD

% [edit] specify some files
file_measurement = "output-measurement_withoutnav-cleaned_2025-02-13-12-47-53.mat";
file_femurSTL    = "Femur_adaptiveremeshed.stl";
file_tibiaSTL    = "Tibia_adaptiveremeshed.stl";

% load the measurements
load(fullfile(path_outputs, file_measurement));
% get the number of groups
n_groups = length(all_measurements);
% [!] Get the number of peaks. Context: There are multiple type of US 
%     measurement point clouds.
%     - data with nav    : cleaned, window-updated
%     - data without nav : cleaned, threshold 1, 2, 3.
n_peaks = size(all_measurements(1).probes(1).peak_2d, 2);

% [edit] Choose the peak type
PEAKIDX_WINAV_CLEAN         = 1;
PEAKIDX_WINAV_WINDOWUPDATE  = 2; 
PEAKIDX_WONAV_USERSELECT    = 1;
PEAKIDX_WONAV_THRESH1       = 2;
PEAKIDX_WONAV_THRESH2       = 3;
PEAKIDX_WONAV_THRESH3       = 4;
peak_selection = PEAKIDX_WONAV_THRESH2;

% check if you select wrong
if(peak_selection > n_peaks)
    error('You selected invalid peak type selection. Check your output-measurement how many types of peak you have.');
end

% create a struct for storing our US point cloud
femurUS = struct('name', "", 'peak_3d', [], 'peak_intensity', [], 'peak_direction', []);
tibiaUS = struct('name', "", 'peak_3d', [], 'peak_intensity', [], 'peak_direction', []);

% obtain the 3d point cloud
for group_idx=1:n_groups
    
    % get the probes measuremnets
    all_probes = all_measurements(group_idx).probes;
    n_probes   = length(all_probes);

    % temporary variables for storing the peak values
    currentgroup_peak3d = [];
    currentgroup_peakintensity = [];
    currentgroup_peakdirection = [];

    % loop for all probe
    for probe_idx=1:n_probes
        % skip if the peak is empty
        if(isempty(all_probes(probe_idx).peak_3d_inref))
            continue;
        end
            
        % get the biggest peak, the first column
        tmp = all_probes(probe_idx).peak_3d_inref(:, peak_selection);
        currentgroup_peak3d = [currentgroup_peak3d, tmp];

        % get the intensity of the biggest peak
        tmp = all_probes(probe_idx).peak_2d(:, peak_selection);
        currentgroup_peakintensity = [currentgroup_peakintensity, tmp];

        % get the "direction" of the points
        tmp1 = all_probes(probe_idx).peak_3d_inref(:, 1);
        tmp2 = all_probes(probe_idx).signal_3d_inref(:, 1);
        tmp3 = (tmp2-tmp1)/norm((tmp2-tmp1));
        currentgroup_peakdirection = [currentgroup_peakdirection, tmp3(1:3)];
    end

    % let's store it into two different group first
    strings = split(all_measurements(group_idx).groupname, '_');
    if (strcmp(strings{2}, 'F'))
        femurUS.name = "femur";
        femurUS.peak_3d = [femurUS.peak_3d, currentgroup_peak3d];
        femurUS.peak_intensity = [femurUS.peak_intensity, currentgroup_peakintensity];
        femurUS.peak_direction = [femurUS.peak_direction, currentgroup_peakdirection];
    elseif(strcmp(strings{2}, 'T'))
        tibiaUS.name = "tibia";
        tibiaUS.peak_3d = [tibiaUS.peak_3d, currentgroup_peak3d];
        tibiaUS.peak_intensity = [tibiaUS.peak_intensity, currentgroup_peakintensity];
        tibiaUS.peak_direction = [tibiaUS.peak_direction, currentgroup_peakdirection];
    end
end

% store them into one variable, easy to access and loopable
allBoneUS  = [femurUS, tibiaUS];
allBoneSTL = {stlread(fullfile(path_bonestl, file_femurSTL)), stlread(fullfile(path_bonestl, file_tibiaSTL))};

% delete some variables
clearvars femurUS tibiaUS ...
          currentgroup_peak3d...
          currentgroup_peakintensity...
          currentgroup_peakdirection...
          tmp tmp1 tmp2 tmp3 strings ...
          probe_idx group_idx all probes;

%% PARAMETERS   

% PCA-alignment parameters
pca_params = struct('name', "", 'correction_R', [], 'correction_t', []);
% femur
pca_params(1).name         = "femur";
pca_params(1).correction_R = eul2rotm([0 pi 0], "ZYX");
pca_params(1).correction_t = [-10 50 -20];
% tibia
pca_params(2).name         = "tibia";
pca_params(2).correction_R = eul2rotm([0 pi 0], "ZYX");
pca_params(2).correction_t = [-20 -10 -40];

% RICP parameters
ricp_params = struct('name', "", 'random_R', 0, 'random_t', 0, 'inlier_ratio', 0, 'icp_iteration', 0, 'ricp_iteration', 0, 'knn_rmse', 0);
% femur
ricp_params(1).name           = "femur";
ricp_params(1).random_R       = 10;
ricp_params(1).random_t       = 10;
ricp_params(1).inlier_ratio   = 0.15;
ricp_params(1).icp_iteration  = 50;
ricp_params(1).ricp_iteration = 100;
ricp_params(1).knn_rmse        = 20;
% tibia
ricp_params(2).name           = "tibia";
ricp_params(2).random_R       = 10;
ricp_params(2).random_t       = 10;
ricp_params(2).inlier_ratio   = 0.20;
ricp_params(2).icp_iteration  = 50;
ricp_params(2).ricp_iteration = 100;
ricp_params(2).knn_rmse       = 20;

% UKF parameters
ukf_params = struct('name', "", 'normal_ratio', 0, 'threshold', 0, 'iteration', 0, 'expected_noise', 0, 'sigma_x_anneal', 0, 'sigma_x_trans', 0, 'sigma_x_theta', 0);
% femur
ukf_params(1).normal_ratio   = 0.5;
ukf_params(1).threshold      = 0.0001;
ukf_params(1).iteration      = 300;
ukf_params(1).expected_noise = 15;
ukf_params(1).sigma_x_anneal = 0.99;
ukf_params(1).sigma_x_trans  = 10;
ukf_params(1).sigma_x_theta  = 10;
% tibia
ukf_params(2).normal_ratio   = 0.5;
ukf_params(2).threshold      = 0.0001;
ukf_params(2).iteration      = 300;
ukf_params(2).expected_noise = 15;
ukf_params(2).sigma_x_anneal = 0.99;
ukf_params(2).sigma_x_trans  = 10;
ukf_params(2).sigma_x_theta  = 10;

% [edit] flag to show bonescans (can't be applied to data without navigation)
is_showbonescan = false;

% [edit] flag to show bone model (original, prereg, ricp, ukf)
is_showpoints   = [false, false, true, true];
is_shownormals  = [false, false, false, false];

% loop for every bone segment
for bone_idx=1:length(allBoneUS)
    %% 1) LOAD DATA

    % get the current boneUS
    boneUS  = allBoneUS(bone_idx);

    % plot the peak3d
    scattersize_all  = boneUS.peak_intensity(2,:);
    scattercolor_all = interp1(linspace(min(scattersize_all), max(scattersize_all), size(colormap1, 1)), colormap1, scattersize_all);
    scatter3( ax1, ...
              boneUS.peak_3d(1,:), ...
              boneUS.peak_3d(2,:), ...
              boneUS.peak_3d(3,:), ...
              30, scattersize_all, 'filled');

    % plot the bone scan
    if(is_showbonescan)

        % get the group names
        groupnames = {all_bonescans.groupname};

        % Extract all group names from the structure and find indices where the 
        % middle part is 'F' for femur and 'T' for tibia
        if(strcmp(boneUS.name, "femur"))
            indices = find(cellfun(@(x) strcmp(x{2}, 'F'), cellfun(@(x) strsplit(x, '_'), groupnames, 'UniformOutput', false)));
        else
            indices = find(cellfun(@(x) strcmp(x{2}, 'T'), cellfun(@(x) strsplit(x, '_'), groupnames, 'UniformOutput', false)));
        end

        % for every group, plot the bone scan
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
    end

    % load bonestl files
    boneCTstl_ref     = allBoneSTL{bone_idx};
    boneCTpc_ref      = pointCloud(boneCTstl_ref.Points);
    boneCTnormals_ref = STLVertexNormals(boneCTstl_ref.ConnectivityList, boneCTstl_ref.Points)';
    boneCTpoints_ref  = [boneCTpc_ref.Location'; ones(1, length(boneCTpc_ref.Location))];
    if(is_showpoints(1))
    scatter3( ax1, ...
              boneCTpoints_ref(1,:), ...
              boneCTpoints_ref(2,:), ...
              boneCTpoints_ref(3,:), ...
              1, "red", "filled", ...
              "Tag", "bone_original");
    end
    if(is_shownormals(1))
    quiver3( ax1, ...
             boneCTpoints_ref(1,:), ...
             boneCTpoints_ref(2,:), ...
             boneCTpoints_ref(3,:), ...
             boneCTnormals_ref(1,:), ...
             boneCTnormals_ref(2,:), ...
             boneCTnormals_ref(3,:), ...
              "Tag", "bone_original");
    end

    % clear unnecessary variables
    clearvars scattersize_all scattercolor_all ...
              groupnames indices group_idx;

    %% 2) PCA ALIGNMENT
    
    % get the principal component of femur stl
    coeff            = pca(boneCTstl_ref.Points);
    mu               = mean(boneCTstl_ref.Points);
    T_pcaBoneCT_ref = [coeff, mu'; 0 0 0 1];
    % display_axis(ax1, T_pcaBoneCT_ref(1:3,4), T_pcaBoneCT_ref(1:3,1:3), 100, "Tag", "baseaxis_boneCT")
    
    % get the principal componen of femur measurement
    coeff           = pca(boneUS.peak_3d(1:3,:)') * pca_params(bone_idx).correction_R;
    mu              = mean(boneUS.peak_3d(1:3,:)') + pca_params(bone_idx).correction_t; 
    T_pcaBoneUS_ref = [coeff, mu'; 0 0 0 1];
    % display_axis(ax1, T_pcaBoneUS_ref(1:3,4), T_pcaBoneUS_ref(1:3,1:3), 100, "Tag", "baseaxis_boneCT");
    
    % transform
    T_pcaBoneUS_pcaBoneCT  = T_pcaBoneUS_ref * inv(T_pcaBoneCT_ref);
    boneCTpoints_pcaBoneCT = T_pcaBoneUS_pcaBoneCT * boneCTpoints_ref;
    % plot
    if(is_showpoints(2))
    scatter3( ax1, ...
              boneCTpoints_pcaBoneCT(1,:), ...
              boneCTpoints_pcaBoneCT(2,:), ...
              boneCTpoints_pcaBoneCT(3,:), ...
              0.5, "magenta", "filled", ...
              "Tag", "bone_pcapoint");
    end
    if(is_shownormals(2))
        % show normals here
    end

    % clear unnecessary variables
    clearvars coeff mu;

    %% 3) R-ICP ALIGNMENT

    pc1 = pointCloud(boneUS.peak_3d(1:3,:)');
    pc2 = pointCloud(boneCTpoints_pcaBoneCT(1:3, :)');
    
    % pc1 and pc2 are defined somewhere above
    iteration = ricp_params(bone_idx).ricp_iteration;
    
    % Preallocate for parfor
    rmse_icp_all = zeros(1, iteration);
    rmse_me_all  = zeros(1, iteration);
    T_icp_all    = cell(1, iteration);
    T_random_all = cell(1, iteration);
    
    for i = 1:iteration
        try
            % make a random transformation
            T_random = rigidtform3d(randomRigidTransfrom(ricp_params(bone_idx).random_R, ricp_params(bone_idx).random_t));
            % apply icp
            [tform, movingReg, rmse_icp] = pcregistericp(pc2, pc1, ...
                "InitialTransform", T_random, ...
                "InlierRatio",      ricp_params(bone_idx).inlier_ratio, ...
                "MaxIterations",    ricp_params(bone_idx).icp_iteration, ...
                "Verbose",          false);
    
            % search nearest neighbor 
            [nn_idx, nn_dist] = knnsearch(movingReg.Location, pc1.Location, "K", ricp_params(bone_idx).knn_rmse);
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
    T_icpBone_pcaBoneUS   = T_icp_all{min_idx};
    % transform
    boneCTpoints_pcaBoneUS  = T_icpBone_pcaBoneUS * T_pcaBoneUS_pcaBoneCT * boneCTpoints_ref;
    boneCTnormals_pcaBoneUS = T_icpBone_pcaBoneUS(1:3, 1:3) * T_pcaBoneUS_pcaBoneCT(1:3, 1:3) * boneCTnormals_ref;
    % plot
    if(is_showpoints(3))
    scatter3( ax1, ...
              boneCTpoints_pcaBoneUS(1,:), ...
              boneCTpoints_pcaBoneUS(2,:), ...
              boneCTpoints_pcaBoneUS(3,:), ...
              0.5, "green", "filled", ...
              "Tag", "bone_ricppoint");
    end
    if(is_shownormals(3))
    quiver3( ax1, ...
             boneCTpoints_pcaBoneUS(1,:), ...
             boneCTpoints_pcaBoneUS(2,:), ...
             boneCTpoints_pcaBoneUS(3,:), ...
             boneCTnormals_pcaBoneUS(1,:), ...
             boneCTnormals_pcaBoneUS(2,:), ...
             boneCTnormals_pcaBoneUS(3,:), ...
              "Tag", "bone_normalpoint");
    end

    % clear unnecessary variables
    clearvars iteration ...
              rmse_icp_all rmse_me_all T_icp_all T_random_all ...
              tform movingReg rmse_icp ...
              nn_idx nn_dist ...
              min_val min_idx;

    %% 4) FINE ALIGNMENT

    [T_icpBone_ukfBone, ~, ~] = ukf_isotropic_registration_ex2( boneUS.peak_3d(1:3,:), boneCTpoints_pcaBoneUS(1:3,:), ...
                                   'movingnormal',  boneUS.peak_direction, ...
                                   'fixednormal',   boneCTnormals_pcaBoneUS, ...
                                   'normalratio',   ukf_params(bone_idx).normal_ratio, ...
                                   'threshold',     ukf_params(bone_idx).threshold, ...
                                   'iteration',     ukf_params(bone_idx).iteration, ...
                                   'expectednoise', ukf_params(bone_idx).expected_noise, ...
                                   'sigmaxanneal',  ukf_params(bone_idx).sigma_x_anneal, ...
                                   'sigmaxtrans',   ukf_params(bone_idx).sigma_x_trans, ...
                                   'sigmaxtheta',   ukf_params(bone_idx).sigma_x_theta, ...
                                   'bestrmse',      true, ...
                                   'verbose',       false, ...
                                   'display',       false);

    % inverse the transformation (because, for some unknown reason, i made the
    % registration from US measurement to the model)
    T_ukfBone_icpBone = inv(T_icpBone_ukfBone);
    % transforrm
    boneCTpoints_final  = T_ukfBone_icpBone * T_icpBone_pcaBoneUS * T_pcaBoneUS_pcaBoneCT * boneCTpoints_ref;
    boneCTnormals_final = T_ukfBone_icpBone(1:3, 1:3) * T_icpBone_pcaBoneUS(1:3, 1:3) * T_pcaBoneUS_pcaBoneCT(1:3, 1:3) * boneCTnormals_ref;
    % plot
    if(is_showpoints(4))
    scatter3( ax1, ...
              boneCTpoints_final(1,:), ...
              boneCTpoints_final(2,:), ...
              boneCTpoints_final(3,:), ...
              0.5, "blue", "filled", ...
              "Tag", "bone_ukfpoint");
    end
    if(is_shownormals(4))
    quiver3( ax1, ...
             boneCTpoints_final(1,:), ...
             boneCTpoints_final(2,:), ...
             boneCTpoints_final(3,:), ...
             boneCTnormals_final(1,:), ...
             boneCTnormals_final(2,:), ...
             boneCTnormals_final(3,:), ...
              "Tag", "bone_ukfnormal");
    end

    % clear unnecessary variables

end








