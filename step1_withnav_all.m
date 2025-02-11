%% PREPARE SOME NECESSARY CONSTANTS

clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_function = "functions";
dir_outputs  = "outputs";
dir_trial    = "trial_0006_Session2_02";

% declare some of the important paths
path_function = fullfile(path_root, dir_function);
path_outputs  = fullfile(path_root, dir_outputs);
path_snapshot = fullfile(path_root, "data", dir_trial, "snapshot");
path_bonescan = fullfile(path_root, "data", dir_trial, "bonescan");
addpath(genpath(path_function));

% declare some of important names
groups = [ "A_F_GTR_reobserved", ...
           "A_F_MID_reobserved", ...
           "A_F_LEP", ...
           "A_F_MEP_reobserved", ...
           "A_T_LEP_reobserved", ...
           "A_T_MEP_reobserved", ...
           "A_T_MID_reobserved", ...
           "A_T_MAL_reobserved"];
n_groups = length(groups);
color  = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};

% create base rotation to trasnform qualisys base vector to matlab
% Qualisys has y direction as up, MATLAB has z direction as up
R_tmp = eul2rotm([0 0 pi/2], "ZYX");
t_tmp = [0 0 0]';
baseRotation_Qualisys2Matlab = [R_tmp, t_tmp; 0 0 0 1];

%% PREPARE THE FIGURES

% get the screen size
scr_size  = get(0, 'ScreenSize');
% prepare the figure object for visualizing 3d signal
fig1 = figure("Name", "3D Signal", "Position", [0 0 scr_size(3)/2, scr_size(4)]);
% create a tab group
tg1 = uitabgroup(fig1);
% create tab and axis
ax1 = {};
for i=1:n_groups
    % make a new tab
    tb1 = uitab(tg1, 'Title', groups(i));
    % set the properties of the axes
    ax_tmp = axes(tb1);
    grid(ax_tmp, "on");
    axis(ax_tmp, "equal");
    title(ax_tmp, "Rigid Bodies");
    view(ax_tmp, [45,15]);
    hold(ax_tmp, "on");
    % store the axes
    ax1{i} = ax_tmp;
end

% prepare the figure object for visualizing 3d signal
fig2 = figure("Name", "2D Signal", "Position", [scr_size(3)/2, 0, scr_size(3)/2, scr_size(4)/2]);
% create a tab group
tg2_parent = uitabgroup(fig2);
% create a tab, each tab will represents each group
tb2_parent = {};
for i=1:n_groups
    % make a new tab
    tb2_parent{i} = uitab(tg2_parent, 'Title', groups(i));
end

% Define custom colormaps
colormap1 = parula(256); % First colormap
colormap2 = bone(256);    % Second colormap

% define the size of scatter dot for 3d signal
scatterdot_scale = 200;

% variable to store the measurements
all_measurements = struct('groupname', "", 'probes', []);
all_bonescans = struct('groupname', "", 'voxels', []);

for group_idx=1:n_groups
    %% READ THE SNAPSHOT DATA (RIGID BODY, WINDOW, AND SIGNAL)
    
    % get the current group
    current_group = groups{group_idx};
    % create the current directory
    directory = fullfile(path_snapshot, current_group);
    % read the rigid body data from the CSV files
    current_grouprigidbody = readCSV_stylusSnapshot(directory);
    % get the number of observation
    n_observation = length(current_grouprigidbody);
    
    % read the window data from the CSV files
    current_groupwindow    = readCSV_windowSnapshot(directory);
    
    % read the ultrasound signal from the current group
    n_probes  = 1;
    n_samples = 3500;
    [current_groupultrasound, timestamps, indexes] = readTIFF_USsignal(directory, n_probes, n_samples);
    
    % preparing data_spec structs which needed by the peak detection algorithm
    data_spec.n_ust     = size(current_groupultrasound, 1);
    data_spec.n_samples = size(current_groupultrasound, 2);
    data_spec.n_frames  = size(current_groupultrasound, 3);
    % preparing us_spec structs which needed by the peak detection algorithm
    us_spec.v_sound     = 1540; % m/s
    us_spec.sample_rate = 50 * 1e6; %Hz
    us_spec.ds          = 1e3 * us_spec.v_sound / (2 * us_spec.sample_rate); % mm
    us_spec.dt          = 1/(us_spec.sample_rate); %s
    us_spec.s_vector    = (1:data_spec.n_samples) .* us_spec.ds; % mm
    us_spec.t_vector    = ((1:data_spec.n_samples) .* us_spec.dt) * 1e6; % mu s
    
    % window selection
    winoffset_mm = 1.5;
    
    %% PREPARE THE FIGURE AND AXES
    
    % create another tab group within the tab
    tg2_child = uitabgroup(tb2_parent{group_idx});
    
    % adjust the axes for figure 2
    ax2 = {};
    for usdata_idx=1:n_observation
        % make a new tab
        tb2_child = uitab(tg2_child, 'Title', ['Shot ' num2str(usdata_idx)]);
        ax_tmp = axes(tb2_child);
    
        % adjust the axis properties
        grid(ax_tmp, 'on');
        grid(ax_tmp, 'minor');
        hold(ax_tmp, 'on');
        xlabel(ax_tmp, 'Depth (mm)');
        ylabel(ax_tmp, 'Amplitude');
        ylim(ax_tmp, [0 5000]);
        title(ax_tmp, ['Snapshot ' num2str(usdata_idx)]);
    
        % store the axes object
        ax2{usdata_idx} = ax_tmp;
    end
    
    %% DISPLAY THE 3D SIGNAL AND 2D SIGNAL

    % Variable to store the probe measurement
    probe = struct('number', 0, 'timestamp', 0, 'signal_2d', [], 'peak_2d', [], 'window_2d', [], 'signal_3d_inref', [], 'peak_3d_inref', []);
    
    % for all the snapshot inside the data struct...
    for rbdata_idx=1:n_observation
    
        % get the trasnsformation of the stylus
        T_stylus_global = current_grouprigidbody(rbdata_idx).Stylus.T;
        R_stylus_global = T_stylus_global(1:3, 1:3);
        t_stylus_global = T_stylus_global(1:3, 4);
    
        % get the transformation of the probe
        T_probe_stylus = current_grouprigidbody(rbdata_idx).Probes(6).T;
        R_probe_stylus = T_probe_stylus(1:3, 1:3);
        t_probe_stylus = T_probe_stylus(1:3, 4);
    
        % get the transformation of the ref (it is in the first data in Probe field)
        T_ref_global = current_grouprigidbody(rbdata_idx).Probes(1).T;    
    
        % get the us signal for this current rigid body and the window
        signal = current_groupultrasound(:,:,rbdata_idx)';
        window = current_groupwindow(rbdata_idx).Window;
    
        % get the signal sample that is inside the window for peak detection
        signal_inwindow_idx = find((us_spec.s_vector > window(2)-winoffset_mm) & (us_spec.s_vector < window(2)+winoffset_mm));
        signal_inwindow_amp = signal(signal_inwindow_idx);
    
        % peak detection
        [peaks_vals, peaks_inwindow_idx] = findpeaks( signal_inwindow_amp, ...
                                                      "SortStr", "descend", ...
                                                      "MinPeakProminence", 40, ...
                                                      "NPeaks", 3);
        % convert peak index in the signal within the window to the actual
        % index in the original signal
        peaks_idx = signal_inwindow_idx(1) + peaks_inwindow_idx;
        peaks_mm  = us_spec.s_vector(peaks_idx);
    
        % remove the near field disturbance
        nearfieldinterference_threshold_mm = 1.5;
        nearfieldinterference_threshold_idx = round(nearfieldinterference_threshold_mm/us_spec.ds);
        signal(1:nearfieldinterference_threshold_idx) = 0;
    
        % downsample the signal
        signal_res       = 0.1; %[mm]
        signal_ticks_mm  = 0:signal_res:us_spec.s_vector(end);
        signal_ticks_idx = knnsearch(us_spec.s_vector', signal_ticks_mm', 'K', 1);
        signal_lowres    = zeros(length(signal_ticks_idx), 1);
        parfor k=1:length(signal_ticks_idx)-1
            signal_lowres(k) = mean(signal(signal_ticks_idx(k):signal_ticks_idx(k+1)));
        end
        
        % get the signal sample that is inside the window
        signal_inwindow_bool = (signal_ticks_mm > window(2)-winoffset_mm) & (signal_ticks_mm < window(2)+winoffset_mm);
        signal_outwindow_idx = find(~signal_inwindow_bool);
        signal_outwindow_mm  = signal_ticks_mm(signal_outwindow_idx);
        signal_inwindow_idx  = find(signal_inwindow_bool);
        signal_inwindow_mm   = signal_ticks_mm(signal_inwindow_idx);
    
        % create the 3d signal
        signal3d_all_inglobal        = [zeros(1, length(signal_lowres)); 
                                        signal_ticks_mm; 
                                        zeros(1, length(signal_lowres)); 
                                        ones(1, length(signal_lowres))];
        signal3d_all_intensity       = round( scatterdot_scale * (signal_lowres / max(signal_lowres)) ) + 1;
    
        % create the 3d signal outside the window
        signal3d_outwindow_inglobal  = [zeros(1, length(signal_outwindow_mm)); 
                                        signal_outwindow_mm; 
                                        zeros(1, length(signal_outwindow_mm)); 
                                        ones(1, length(signal_outwindow_mm))];
        signal3d_outwindow_intensity = signal3d_all_intensity(signal_outwindow_idx);
    
        % create the 3d signal within the window
        signal3d_inwindow_inglobal   = [zeros(1, length(signal_inwindow_mm)); 
                                        signal_inwindow_mm; 
                                        zeros(1, length(signal_inwindow_mm)); 
                                        ones(1, length(signal_inwindow_mm))];
        signal3d_inwindow_intensity  = signal3d_all_intensity(signal_inwindow_idx);
    
        % create the 3d signal peak
        peak3d_inglobal              = [zeros(1, length(peaks_mm));
                                        peaks_mm;
                                        zeros(1, length(peaks_mm));
                                        ones(1, length(peaks_mm))];
        
        % transform the signal
        signal3d_all_inref         = inv(T_ref_global) * T_stylus_global * T_probe_stylus * baseRotation_Qualisys2Matlab * signal3d_all_inglobal;
        signal3d_outwindow_inref   = inv(T_ref_global) * T_stylus_global * T_probe_stylus * baseRotation_Qualisys2Matlab * signal3d_outwindow_inglobal;
        signal3d_inwindow_inref    = inv(T_ref_global) * T_stylus_global * T_probe_stylus * baseRotation_Qualisys2Matlab * signal3d_inwindow_inglobal;
        peak3d_inref               = inv(T_ref_global) * T_stylus_global * T_probe_stylus * baseRotation_Qualisys2Matlab * peak3d_inglobal;
    
        % plot the 3d signal
        size_outWindow  = signal3d_outwindow_intensity;
        color_outWindow = interp1(linspace(min(size_outWindow), max(size_outWindow), size(colormap1, 1)), colormap1, size_outWindow);
        scatter3( ax1{group_idx}, ...
                  signal3d_outwindow_inref(1,:), ...
                  signal3d_outwindow_inref(2,:),  ...
                  signal3d_outwindow_inref(3,:), ...
                  size_outWindow, color_outWindow, ...
                  'filled', 'MarkerEdgeColor', 'black', ...
                  'MarkerFaceAlpha', .1, ...
                  'MarkerEdgeAlpha', .2, ...
                  'Tag', "amode_3d");
        size_inWindow  = signal3d_inwindow_intensity;
        color_inWindow = 'red';
        scatter3( ax1{group_idx}, ...
                  signal3d_inwindow_inref(1,:), ...
                  signal3d_inwindow_inref(2,:),  ...
                  signal3d_inwindow_inref(3,:), ...
                  size_inWindow, color_inWindow, ...
                  'filled', 'MarkerEdgeColor', 'red', ...
                  'MarkerFaceAlpha', .15, ...
                  'MarkerEdgeAlpha', .15, ...
                  'Tag', "amode_3d");
        if(~isempty(peak3d_inref))
        scatter3( ax1{group_idx}, ...
                  peak3d_inref(1), ...
                  peak3d_inref(2), ...
                  peak3d_inref(3), ...
                  50, 'yellow', ...
                  'filled', 'MarkerEdgeColor', 'black', ...
                  'LineWidth', 1.5, ...
                  'Tag', "amode_3d");
        end
        text( ax1{group_idx},...
              signal3d_outwindow_inref(1,1), ...
              signal3d_outwindow_inref(2,1), ...
              signal3d_outwindow_inref(3,1), ...
              string(rbdata_idx), ...
              "HorizontalAlignment", "center", ...
              "VerticalAlignment", "middle");
    
        % plot the 2d signal
        plot(ax2{rbdata_idx}, us_spec.s_vector, current_groupultrasound(1,:,rbdata_idx), "Color", 'b', "LineWidth", 1);
        if(~isempty(peaks_mm)) plot(ax2{rbdata_idx}, peaks_mm, peaks_vals, "or", "MarkerFaceColor", "r"); end
        xline(ax2{rbdata_idx}, window(2)-winoffset_mm, "Color", 'r', "LineWidth", 2);
        xline(ax2{rbdata_idx}, window(2)+winoffset_mm, "Color", 'r', "LineWidth", 2);

        probe(rbdata_idx).number    = rbdata_idx;
        probe(rbdata_idx).timestamp = current_grouprigidbody.Timestamp;
        probe(rbdata_idx).signal_2d = [us_spec.s_vector; signal'];
        probe(rbdata_idx).peak_2d   = [peaks_mm; peaks_vals'];
        probe(rbdata_idx).window_2d = [window(2)-winoffset_mm window(2) window(2)+winoffset_mm];
        probe(rbdata_idx).signal_3d_inref = signal3d_all_inref;
        probe(rbdata_idx).peak_3d_inref   = peak3d_inref;

    end
    
    %% DISPLAY THE SCANNED BONE SURFACE
    
    tmp = strsplit(current_group, '_');
    current_group = strjoin(tmp(1:3), '_');
    
    % load the lookup table for surfaces
    filename_surfacecsv = "surfaces_processed.csv";
    fullfile_surfacecsv = fullfile(path_bonescan, filename_surfacecsv);
    lookuptbl_surfacecsv = readtable(fullfile_surfacecsv, "Delimiter", ',');
    
    % find the rows where the 'Area' column matches what we are interested in
    areaColumn   = string(lookuptbl_surfacecsv.Area);
    matchingRows = lookuptbl_surfacecsv(areaColumn == current_group, :);
    filename_surface = string(matchingRows.Name);
    fullfile_surface = fullfile(path_bonescan, filename_surface);
    
    % read the surface
    reader = MHAReader(fullfile_surface);
    if ~reader.readVolumeImage()
        error('FileOpenError:CannotOpen', 'Failed to open the file: %s', fullfile_surface);
    end
    headerInfo = reader.getMHAHeader();
    volumeData = reader.getMHAVolume();
    
    % find indices over threshold
    threshold = matchingRows.Threshold;
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
    
    % Normalize z values to range [0, 1] for colormap mapping
    z = voxelcoordinate_homogeneous(3,:);
    colors2 = interp1(linspace(min(z), max(z), size(colormap2, 1)), colormap2, z);
    
    % display the points
    scatter3(  ax1{group_idx}, ...
               voxelcoordinate_homogeneous(1,:), ...
               voxelcoordinate_homogeneous(2,:),  ...
               voxelcoordinate_homogeneous(3,:), ...
               5, colors2, ...
              'filled', ...
              'Tag', "bone_surface");

    all_measurements(group_idx).groupname = current_group;
    all_measurements(group_idx).probes    = probe;

    all_bonescans(group_idx).groupname = current_group;
    all_bonescans(group_idx).voxels = voxelcoordinate_homogeneous;
end

% currentTime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% file_output = ['measurementoutput_' currentTime '.mat'];
% save(fullfile(path_outputs, file_output), 'all_measurements', 'all_bonescans');










