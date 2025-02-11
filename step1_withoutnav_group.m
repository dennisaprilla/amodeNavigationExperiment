%% PREPARE SOME NECESSARY CONSTANTS

clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_function = "functions";
dir_trial    = "trial_0003_Session1_02";

% declare some of the important paths
path_function = fullfile(path_root, dir_function);
path_snapshot = fullfile(path_root, "data", dir_trial, "snapshot");
path_bonescan = fullfile(path_root, "data", dir_trial, "bonescan");
addpath(genpath(path_function));

% declare some of important names
groups = [ "A_F_GTR", ...
           "A_F_MID", ...
           "A_F_LEP", ...
           "A_F_MEP", ...
           "A_T_LEP", ...
           "A_T_MEP", ...
           "A_T_MID", ...
           "A_T_MAL"];
color  = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};

% create base rotation to trasnform qualisys base vector to matlab
% Qualisys has y direction as up, MATLAB has z direction as up
R_tmp = eul2rotm([0 0 pi/2], "ZYX");
t_tmp = [0 0 0]';
baseRotation_Qualisys2Matlab = [R_tmp, t_tmp; 0 0 0 1];

%% READ THE SNAPSHOT DATA (RIGID BODY, WINDOW, AND SIGNAL)

% select the group
group_idx = 1;
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

% normalize ultrasound (?) i don't know that this is necessary but

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
winoffset_mm = 1.10;

%% PREPARE THE FIGURE AND AXES

% get the screen size
scr_size  = get(0, 'ScreenSize');
% prepare the figure object
fig1 = figure("Name", "3D Signal", "Position", [0 0 scr_size(3)/2, scr_size(4)]);
fig2 = figure("Name", "2D Signal", "Position", [scr_size(3)/2, 0, scr_size(3)/2, scr_size(4)/2]);
% for fig2, i will need to show the tabs for amode 2d signals
tab_group = uitabgroup(fig2);

% adjust the axis for figure 1
ax1  = axes(fig1);
grid(ax1, "on");
axis(ax1, "equal");
title(ax1, "Rigid Bodies");
view(ax1, [45,15]);
hold(ax1, "on");

% adjust the axes for figure 2
ax2 = {};
for usdata_idx=1:n_observation
    % make a new tab
    tab = uitab(tab_group, 'Title', ['Shot ' num2str(usdata_idx)]);
    ax_tmp = axes(tab);

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

% Define custom colormaps
colormap1 = parula(256); % First colormap
colormap2 = bone(256);    % Second colormap

% define the size of scatter dot for 3d signal
scatterdot_scale = 200;

%% DISPLAY THE 3D SIGNAL AND 2D SIGNAL

% for all the snapshot inside the data struct...
for rbdata_idx=1:n_observation
    %% 1) GET ALL THE NECESSARY DATA

    % 1.1) get the trasnsformation of the stylus
    T_stylus_global = current_grouprigidbody(rbdata_idx).Stylus.T;
    R_stylus_global = T_stylus_global(1:3, 1:3);
    t_stylus_global = T_stylus_global(1:3, 4);

    % 1.2) get the transformation of the probe
    T_probe_stylus = current_grouprigidbody(rbdata_idx).Probes(6).T;
    R_probe_stylus = T_probe_stylus(1:3, 1:3);
    t_probe_stylus = T_probe_stylus(1:3, 4);

    % 1.3) get the transformation of the ref (it is in the first data in Probe field)
    T_ref_global = current_grouprigidbody(rbdata_idx).Probes(1).T;    

    % 1.4) get the us signal for this current rigid body and the window
    signal = current_groupultrasound(:,:,rbdata_idx)';
    window = current_groupwindow(rbdata_idx).Window;

    %% 2) PEAK DETECTION FROM WINDOW

    % 2.1) get the signal sample that is inside the window for peak detection
    signal_inwindow_idx = find((us_spec.s_vector > window(2)-winoffset_mm) & (us_spec.s_vector < window(2)+winoffset_mm));
    signal_inwindow_amp = signal(signal_inwindow_idx);

    % 2.2) peak detection
    [peak_fromwindow_val, peak_inwindow_idx] = findpeaks( signal_inwindow_amp, ...
                                               "SortStr", "descend", ...
                                               "NPeaks", 1);
    % 2.3) convert peak index in the signal within the window to the actual
    %      index in the original signal
    peak_fromwindow_idx = signal_inwindow_idx(1) + peak_inwindow_idx;
    peak_fromwindow_mm  = us_spec.s_vector(peak_fromwindow_idx);

    %% 3) PEAK DETECTION FROM (MULTIPLE) THRESHOLDS

    % determine the noise
    mean_noise       = mean(signal(end-100:end));
    thresh_constants = [10, 20, 30]; 

    % variable to stores the peak
    peak_fromthresh_idx = [];
    peak_fromthresh_mm  = [];
    peak_fromthresh_amp = [];
    % peak detection with multiple thresholds
    for thresh_constant=thresh_constants
        [tmp_amp, tmp_idcs] = findpeaks( signal_inwindow_amp, ...
                                   "SortStr", "none", ...
                                   "MinPeakHeight", thresh_constant*mean_noise, ...
                                   "MinPeakProminence", 5*mean_noise);
        tmp_idx = tmp_idcs(end);
        tmp_mm = us_spec.s_vector(peak_fromthresh1_idx);

        peak_fromthresh_idx = [peak_fromthresh_idx, tmp_idx];
        peak_fromthresh_mm  = [peak_fromthresh_mm, tmp_mm];
        peak_fromthresh_amp = [peak_fromthresh_amp, tmp_amp];
    end

    %% 4) MAKING 3D SIGNAL
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
    peak3d_fromwindow_inglobal   = [0, peak_fromwindow_mm, 0, 1]';
    peak3d_fromthresh_inglobal   = [zeros(1, length(peak_fromthresh_mm));
                                    peak_fromthresh1_mm;
                                    zeros(1, length(peak_fromthresh_mm));
                                    ones(1, length(peak_fromthresh_mm))];
    
    % transform the signal
    T_all = inv(T_ref_global) * T_stylus_global * T_probe_stylus * baseRotation_Qualisys2Matlab;
    signal3d_outwindow_inprobe   = T_all * signal3d_outwindow_inglobal;
    signal3d_inwindow_inprobe    = T_all * signal3d_inwindow_inglobal;
    peak3d_fromwindow_inprobe    = T_all * peak3d_fromwindow_inglobal;
    peak3d_fromthresh_inprobe    = T_all * peak3d_fromthresh_inglobal;

    % plot the 3d signal
    size_outWindow  = signal3d_outwindow_intensity;
    color_outWindow = interp1(linspace(min(size_outWindow), max(size_outWindow), size(colormap1, 1)), colormap1, size_outWindow);
    scatter3( ax1, ...
              signal3d_outwindow_inprobe(1,:), ...
              signal3d_outwindow_inprobe(2,:),  ...
              signal3d_outwindow_inprobe(3,:), ...
              size_outWindow, color_outWindow, ...
              'filled', 'MarkerEdgeColor', 'black', ...
              'MarkerFaceAlpha', .1, ...
              'MarkerEdgeAlpha', .2, ...
              'Tag', "amode_3d");
    size_inWindow  = signal3d_inwindow_intensity;
    color_inWindow = 'red';
    scatter3( ax1, ...
              signal3d_inwindow_inprobe(1,:), ...
              signal3d_inwindow_inprobe(2,:),  ...
              signal3d_inwindow_inprobe(3,:), ...
              size_inWindow, color_inWindow, ...
              'filled', 'MarkerEdgeColor', 'red', ...
              'MarkerFaceAlpha', .15, ...
              'MarkerEdgeAlpha', .15, ...
              'Tag', "amode_3d");
    scatter3( ax1, ...
              peak3d_fromwindow_inprobe(1), ...
              peak3d_fromwindow_inprobe(2), ...
              peak3d_fromwindow_inprobe(3), ...
              50, 'yellow', ...
              'filled', 'MarkerEdgeColor', 'black', ...
              'LineWidth', 1.5, ...
              'Tag', "amode_3d");
    text( ax1,...
          signal3d_outwindow_inprobe(1,1), ...
          signal3d_outwindow_inprobe(2,1), ...
          signal3d_outwindow_inprobe(3,1), ...
          string(rbdata_idx), ...
          "HorizontalAlignment", "center", ...
          "VerticalAlignment", "middle");

    % plot the 2d signal
    plot(ax2{rbdata_idx}, us_spec.s_vector, current_groupultrasound(1,:,rbdata_idx), "Color", 'b', "LineWidth", 1);
    plot(ax2{rbdata_idx}, peak_fromwindow_mm, peak_val, "or", "MarkerFaceColor", "r");
    xline(ax2{rbdata_idx}, window(2)-winoffset_mm, "Color", 'r', "LineWidth", 2);
    xline(ax2{rbdata_idx}, window(2)+winoffset_mm, "Color", 'r', "LineWidth", 2);
end

%% DISPLAY THE SCANNED BONE SURFACE

% this script is intended for quantifying the quality of the measurement
% without navigation system








