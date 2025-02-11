clear; clc; close all;

directory = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a\data\trial_0006_Session2_02\snapshot\A_F_GTR";

% Get a list of all CSV files in the folder
csvFiles = dir(fullfile(directory, '*_WindowConfig.csv'));

% Preallocate a structure array
allDataStruct = struct('Timestamp', {}, 'Window', []);

% Loop through each file and read it
for i = 1:length(csvFiles)
    % Extract the filename
    filename = csvFiles(i).name;

    % Extract the timestamp from the filename
    % Assuming the timestamp is the first part before "_"
    timestampStr = regexp(filename, '^\d+', 'match', 'once');
    if isempty(timestampStr)
        warning('Filename "%s" does not start with a valid timestamp. Skipping file.', filename);
        continue;
    end
    timestamp = str2double(timestampStr);

    % Generate the full file path
    filePath = fullfile(csvFiles(i).folder, csvFiles(i).name);

    % Read the CSV file into a table with error handling
    try
        tempData = readtable(filePath);
    catch ME
        warning('Failed to read file "%s": %s. Skipping file.', filename, ME.message);
        continue;
    end

    % Extract the row with "B_N_STY" for its GroupName column
    targetRows = strcmp(tempData.GroupName, 'B_N_STY');
    stylusRow  = tempData(targetRows, :);

    % Store timestamp and window
    allDataStruct(i).Timestamp = timestamp;
    allDataStruct(i).Window = [stylusRow.LowerBound, stylusRow.Middle, stylusRow.UpperBound];

    % Optionally, display the name of the file read
    fprintf('Processed file: %s with timestamp %d\n', filename, timestamp);
end