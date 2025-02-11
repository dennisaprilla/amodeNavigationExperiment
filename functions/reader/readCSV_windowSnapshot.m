function allDataStruct = readCSV_windowSnapshot(directory)
% READCSV_WINDOWSNAPSHOT Reads and processes CSV files containing window configuration data.
%
% This function processes all CSV files in the specified directory whose 
% names end with '_WindowConfig.csv'. It extracts timestamp information 
% from the filenames, reads the data from the CSV files, and collects 
% specific rows where the 'GroupName' column matches 'B_N_STY'. For each 
% matching row, it stores the timestamp and the window configuration 
% (LowerBound, Middle, UpperBound) into a structured array.
%
% INPUT:
%   directory (string): Path to the directory containing the CSV files.
%
% OUTPUT:
%   allDataStruct (struct array): A structure array containing the following fields:
%       - Timestamp: The timestamp extracted from the filename.
%       - Window: A vector containing [LowerBound, Middle, UpperBound] from the row
%                 where 'GroupName' is 'B_N_STY'.
%
% The function handles errors gracefully by skipping files that do not
% conform to the expected format or are not readable, and issues warnings
% when such issues are encountered.
%
% Notes:
% - Filenames are expected to start with a numeric timestamp followed by '_WindowConfig.csv'.
% - If a file does not contain a 'B_N_STY' row in the 'GroupName' column, 
%   it will still be processed, but an empty Window will be stored.
%
% Example:
%   result = readCSV_windowSnapshot('C:\Data\WindowConfigs');
%
%   This will read all '_WindowConfig.csv' files in the directory
%   'C:\Data\WindowConfigs' and return the structured array `result` containing
%   timestamps and window data.

% Get a list of all CSV files in the folder
csvFiles = dir(fullfile(directory, '*_WindowConfig.csv'));

% Preallocate a structure array
allDataStruct = struct('Timestamp', {}, 'Window', []);

% extra offset for window
offset_mm = 0;

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
    allDataStruct(i).Window = [stylusRow.LowerBound-offset_mm, stylusRow.Middle, stylusRow.UpperBound+offset_mm];

    % Optionally, display the name of the file read
    fprintf('Processed file: %s with timestamp %d\n', filename, timestamp);
end

end

