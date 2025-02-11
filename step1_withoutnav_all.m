%% PREPARE SOME NECESSARY CONSTANTS

clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_function = "functions";
dir_outputs  = "outputs";
dir_trial    = "trial_0003_Session1_02";

% declare some of the important paths
path_function = fullfile(path_root, dir_function);
path_outputs  = fullfile(path_root, dir_outputs);
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
n_groups = length(groups);
% Threshold constants
thresh_allconstants = [ 20, 10, 5;
                        12,  7, 3;
                        15, 10, 5;
                        30, 20, 10;
                        30, 20, 10;
                        30, 20, 10;
                        30, 20, 10;
                        30, 20, 10];

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

%% THIS IS JUST ADDITIONAL
% In some measurements, i noticed that T_ref_global is not being tracked by
% the motion capture system. It will gives NaN. Of course it will be
% disasterous for the algorithms. So, since this particular experiment
% assumes that T_ref_global is static (moves only because of noise), so i
% will just average (and also properly orthogonalize again) all the
% detected T_ref_global.

% variable that stores T_ref_global, for emergency if there are
% measurement that does not detect the T_ref_global
T_ref_global_tmp = eye(4);

% Variable to store R_tmp and t_tmp, for emergency if there are
% measurement that does not detect the T_ref_global
R_sum = zeros(3,3);
t_sum = zeros(3,1);
count = 0;

fprintf("Averaging T_ref_global...\n");
for group_idx=1:n_groups
    % get the current group
    current_group = groups{group_idx};
    % create the current directory
    directory = fullfile(path_snapshot, current_group);
    % read the rigid body data from the CSV files
    current_grouprigidbody = readCSV_stylusSnapshot(directory);
    % get the number of observation
    n_observation = length(current_grouprigidbody);

    % loop through all observations
    for rbdata_idx=1:n_observation
        % get the transformation of the ref (it is in the first data in Probe field)
        T_ref_global = current_grouprigidbody(rbdata_idx).Probes(1).T;
        
        % This part is extra, i just want to store and average the
        % T_ref_global, in case there are some measurement that missing it
        if(~any(isnan(T_ref_global), 'all'))
            R_sum = R_sum+T_ref_global(1:3, 1:3);
            t_sum = t_sum+T_ref_global(1:3, 4);
            count = count+1;
        end
    end
end

% Compute the arithmetic means
R_mean = R_sum / count;
t_avg = t_sum / count;

% Project R_mean to the closest rotation matrix via SVD
[U, ~, V] = svd(R_mean);
R_avg = U * V';

% Assemble the average transformation matrix
T_ref_global_avg = eye(4);
T_ref_global_avg(1:3, 1:3) = R_avg;
T_ref_global_avg(1:3, 4) = t_avg;

fprintf("Averaging T_ref_global done.");
pause(2);
clc;

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
        if(any(isnan(T_ref_global), 'all'))
            T_ref_global = T_ref_global_avg;
        end
    
        % 1.4) get the us signal for this current rigid body and the window
        signal = current_groupultrasound(:,:,rbdata_idx)';
        window = current_groupwindow(rbdata_idx).Window;

        % remove the near field disturbance
        nearfieldinterference_threshold_mm = 1.5;
        nearfieldinterference_threshold_idx = round(nearfieldinterference_threshold_mm/us_spec.ds);
        signal(1:nearfieldinterference_threshold_idx) = 0;

        %% 2) PEAK DETECTION FROM WINDOW
    
        % 2.1) get the signal sample that is inside the window for peak detection
        signal_inwindow_idx = find((us_spec.s_vector > window(2)-winoffset_mm) & (us_spec.s_vector < window(2)+winoffset_mm));
        signal_inwindow_amp = signal(signal_inwindow_idx);
    
        % 2.2) peak detection
        [peaks_fromwindow_vals, peaks_inwindow_idx] = findpeaks( signal_inwindow_amp, ...
                                                      "SortStr", "descend", ...
                                                      "MinPeakProminence", 40, ...
                                                      "NPeaks", 1);
        % 2.2) convert peak index in the signal within the window to the actual
        %      index in the original signal
        peaks_fromwindow_idx = signal_inwindow_idx(1) + peaks_inwindow_idx;
        peaks_fromwindow_mm  = us_spec.s_vector(peaks_fromwindow_idx);

        %% 3) PEAK DETECTION FROM (MULTIPLE) THRESHOLDS

        % determine the noise
        mean_noise       = mean(signal(end-100:end));
        thresh_constants = thresh_allconstants(group_idx,:); 
    
        % variable to stores the peak
        peaks_fromthresh_idx = [];
        peaks_fromthresh_mm  = [];
        peaks_fromthresh_amp = [];

        % peak detection with multiple thresholds
        for thresh_constant=thresh_constants
            [peak_amps, peak_idcs] = findpeaks( signal, ...
                                             "SortStr", "none", ...
                                             "MinPeakHeight", thresh_constant*mean_noise, ...
                                             "MinPeakProminence", thresh_constant*mean_noise, ...
                                             "MinPeakDistance", round(1.5/us_spec.ds)); %2mm

            % Check the adjacent peak (left), if that peak is too close a
            thresh_dist  = 2; %mm
            thresh_ratio = 1.75;

            % let's check every peak, but let's begin with selecting the
            % last peak as our focus peak.
            sel_idx_peak_idcs = length(peak_idcs);
            prev_sel_idx_peak_idcs = NaN;
            while(sel_idx_peak_idcs>1)
                % get the difference between every peak to the last (selected) peak
                current_mm  = us_spec.s_vector(peak_idcs(sel_idx_peak_idcs));
                allleft_mm  = us_spec.s_vector(peak_idcs(1:sel_idx_peak_idcs-1));
                diffs_mm    = abs(allleft_mm-current_mm);
                % from the diffs, get the indices that is below the threshold distance
                belowdistthresh_idcs_peak_idcs = find(diffs_mm<thresh_dist);
                % if there is none, we break the loop, means that we select
                % this particular peak
                if(isempty(belowdistthresh_idcs_peak_idcs))
                    break;
                end                

                % now, probably there are more than one peaks that falls
                % within distance threshold, we should check one by one,
                % starting from the right most one
                pot_idx_belowthresh_idcs = length(belowdistthresh_idcs_peak_idcs);
                while(pot_idx_belowthresh_idcs>0)
                    % get the potential selected peak
                    pot_idx_peak_idcs = belowdistthresh_idcs_peak_idcs(pot_idx_belowthresh_idcs);
                    % calculate the ratio between the potential selected peak
                    % amplitude with the
                    current_ratio = signal(peak_idcs(pot_idx_peak_idcs)) / signal(peak_idcs(sel_idx_peak_idcs));
                    % if the ratio is higher than the selected peak, the
                    % potential peak will be the new peak
                    if(current_ratio>thresh_ratio)
                        sel_idx_peak_idcs = pot_idx_peak_idcs;
                    end
                    
                    % go to the next peak on the left
                    pot_idx_belowthresh_idcs = pot_idx_belowthresh_idcs-1;
                end

                % we should break the loop if we don't find any peak that
                % is below the distance threshold with ratio bigger than
                % ratio treshold. This can be known by the previous
                % selected peak in the loop is the same as the current one
                if(prev_sel_idx_peak_idcs==sel_idx_peak_idcs)
                    break;
                % if not, we still need to check if there is other
                else
                    prev_sel_idx_peak_idcs = sel_idx_peak_idcs;
                end               
                
            end
            
            % get the properties of the current peak 
            peak_idx = peak_idcs(sel_idx_peak_idcs);
            peak_mm  = us_spec.s_vector(peak_idx);
            peak_amp = signal(peak_idx);
    
            % store it to the respective variables
            peaks_fromthresh_idx = [peaks_fromthresh_idx, peak_idx];
            peaks_fromthresh_mm  = [peaks_fromthresh_mm, peak_mm];
            peaks_fromthresh_amp = [peaks_fromthresh_amp, peak_amp];
        end

        %% 4) MAKING 3D SIGNAL
    
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
        peak3d_fromwindow_inglobal   = [zeros(1, length(peaks_fromwindow_mm));
                                        peaks_fromwindow_mm;
                                        zeros(1, length(peaks_fromwindow_mm));
                                        ones(1, length(peaks_fromwindow_mm))];
        peak3d_fromthresh_inglobal   = [zeros(1, length(peaks_fromthresh_mm));
                                        peaks_fromthresh_mm;
                                        zeros(1, length(peaks_fromthresh_mm));
                                        ones(1, length(peaks_fromthresh_mm))];
        
        % transform the signal
        T_all = inv(T_ref_global) * T_stylus_global * T_probe_stylus * baseRotation_Qualisys2Matlab;
        signal3d_all_inref         = T_all * signal3d_all_inglobal;
        signal3d_outwindow_inref   = T_all * signal3d_outwindow_inglobal;
        signal3d_inwindow_inref    = T_all * signal3d_inwindow_inglobal;
        peak3d_fromwindow_inref    = T_all * peak3d_fromwindow_inglobal;
        peak3d_fromthresh_inref    = T_all * peak3d_fromthresh_inglobal;
    
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
        scatter3( ax1{group_idx}, ...
                  peak3d_fromthresh_inref(1,:), ...
                  peak3d_fromthresh_inref(2,:), ...
                  peak3d_fromthresh_inref(3,:), ...
                  40, 'blue', 'filled', ...
                  'MarkerEdgeColor', 'black', ...
                  'LineWidth', 1, ...
                  'Tag', "amode_3d");
        if(~isempty(peak3d_fromwindow_inref))
        scatter3( ax1{group_idx}, ...
                  peak3d_fromwindow_inref(1), ...
                  peak3d_fromwindow_inref(2), ...
                  peak3d_fromwindow_inref(3), ...
                  100, "yellow", 'filled', ...
                  'MarkerEdgeColor', 'black', ...
                  'LineWidth', 1, ...
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
        if(~isempty(peaks_fromwindow_mm)) 
            plot(ax2{rbdata_idx}, peaks_fromwindow_mm, peaks_fromwindow_vals, "ob", "MarkerFaceColor", "y", "MarkerSize", 12); 
        end
        plot(ax2{rbdata_idx}, peaks_fromthresh_mm, peaks_fromthresh_amp, "ob", "MarkerFaceColor", "b"); 
        xline(ax2{rbdata_idx}, window(2)-winoffset_mm, "Color", 'r', "LineWidth", 2);
        xline(ax2{rbdata_idx}, window(2)+winoffset_mm, "Color", 'r', "LineWidth", 2);

        probe(rbdata_idx).number          = rbdata_idx;
        probe(rbdata_idx).timestamp       = current_grouprigidbody.Timestamp;
        probe(rbdata_idx).signal_2d       = [us_spec.s_vector; signal'];
        probe(rbdata_idx).peak_2d         = [peaks_fromwindow_mm, peaks_fromthresh_mm; peaks_fromwindow_vals, peaks_fromthresh_amp];
        probe(rbdata_idx).window_2d       = [window(2)-winoffset_mm window(2) window(2)+winoffset_mm];
        probe(rbdata_idx).signal_3d_inref = signal3d_all_inref;
        probe(rbdata_idx).peak_3d_inref   = [peak3d_fromwindow_inref, peak3d_fromthresh_inref];

    end

    %% DISPLAY THE SCANNED BONE SURFACE
    
    % this script is intended for quantifying the quality of the measurement
    % without navigation system

    all_measurements(group_idx).groupname = current_group;
    all_measurements(group_idx).probes    = probe;

    all_bonescans(group_idx).groupname = current_group;
    all_bonescans(group_idx).voxels = [];

end

currentTime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
file_output = ['measurementoutput_withoutnav_' currentTime '.mat'];
save(fullfile(path_outputs, file_output), 'all_measurements', 'all_bonescans');










