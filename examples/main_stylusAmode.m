clear; clc; close all;

% change this
path_root    = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a";
dir_function = "functions";
dir_trial    = "trial_0006_Session2_02";

% declare some of the important paths
path_function = fullfile(path_root, dir_function);
path_snapshot = fullfile(path_root, "data", dir_trial, "snapshot");
addpath(genpath(path_function));

% declare some of important names
groups = {"A_F_GTR", "A_F_MID", "A_F_LEP", "A_F_MEP", "A_T_LEP", "A_T_MEP", "A_T_MID", "A_T_MAL"};

% for each group...
for group_idx=1:length(groups)
    % get the current group
    current_group = groups{group_idx};
    % create the current directory
    directory = fullfile(path_snapshot, current_group);

    % read the ultrasound signal from the current group
    n_probes  = 1;
    n_samples = 3500;
    [all_ultrasoundfrd, timestamps, indexes] = readTIFF_USsignal(directory, n_probes, n_samples);
    
    % preparing data_spec structs which needed by the peak detection algorithm
    data_spec.n_ust     = size(all_ultrasoundfrd, 1);
    data_spec.n_samples = size(all_ultrasoundfrd, 2);
    data_spec.n_frames  = size(all_ultrasoundfrd, 3);
    % preparing us_spec structs which needed by the peak detection algorithm
    us_spec.v_sound     = 1540; % m/s
    us_spec.sample_rate = 50 * 1e6; %Hz
    us_spec.ds          = 1e3 * us_spec.v_sound / (2 * us_spec.sample_rate); % mm
    us_spec.dt          = 1/(us_spec.sample_rate); %s
    us_spec.d_vector    = (1:data_spec.n_samples) .* us_spec.ds; % mm
    us_spec.t_vector    = ((1:data_spec.n_samples) .* us_spec.dt) * 1e6; % mu s
    
    % Create a figure with tabs to select measurement
    f = figure;
    tab_group = uitabgroup(f);
    n_measurements = size(all_ultrasoundfrd, 3);
    
    % Create tabs dynamically for each measurement
    for i = 1:n_measurements
        % make a new tab
        tab = uitab(tab_group, 'Title', ['Shot ' num2str(i)]);
        ax = axes(tab);
        
        % Plot the corresponding measurement
        plot(ax, us_spec.d_vector, all_ultrasoundfrd(1,:,i));
        grid(ax, 'on');
        grid(ax, 'minor');
        xlabel(ax, 'Depth (mm)');
        ylabel(ax, 'Amplitude');
        title(ax, ['Snapshot ' num2str(i)]);
    end
    
end



