function allDataStruct = readCSV_stylusSnapshot(directory)
%readCSV_stylusSnapshot Reads stylus and probes data from CSV files in a directory.
%
%   allDataStruct = readCSV_stylusSnapshot(directory) reads all CSV files in the 
%   specified directory, extracts timestamp, stylus, and probes data, and returns 
%   a structure array containing this information for each file.
%
%   Input:
%       directory - A string specifying the path to the folder containing the CSV files.
%                  The directory should contain files with the expected format.
%
%   Output:
%       allDataStruct - A structure array where each element corresponds to one CSV file.
%                       Each element contains the following fields:
%                       
%           Timestamp : Numeric value representing the timestamp extracted from the filename.
%           Stylus    : A struct with the following fields:
%                       - name: The name of the stylus (string).
%                       - T   : A 4x4 transformation matrix representing the stylus's position and orientation.
%           Probes    : An array of structs, where each struct has the following fields:
%                       - name: The name of the probe (string).
%                       - T   : A 4x4 transformation matrix representing the probe's position and orientation.
%
% USAGE:
%   allDataStruct = readCSV_stylusSnapshot('path/to/file.csv');
%
% EXAMPLE:
%   allDataStruct = readCSV_stylusSnapshot('MyData.csv');
%   disp(allDataStruct);

% Get a list of all CSV files in the folder
csvFiles = dir(fullfile(directory, '*_MocapRecording.csv'));

% Preallocate a structure array
allDataStruct = struct('Timestamp', {}, 'Stylus', {}, 'Probes', {});

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
    
    % Extract the first row for Stylus data
    stylusRow = tempData(1, :);

    % Matlab string can be weird, better check the datatype first
    if iscell(stylusRow.name)
        stylusName = stylusRow.name{1};
    elseif ischar(stylusRow.name) || isstring(stylusRow.name)
        stylusName = stylusRow.name;
    else
        error('Unexpected data type for "name" in file: %s', filename);
    end

    % Get the quaternion and the translation
    qStylus = [stylusRow.q4, stylusRow.q1, stylusRow.q2, stylusRow.q3];     % Quaternion (q4 is the scalar, i put it in the front using matlab notation)
    tStylus = [stylusRow.t1, stylusRow.t2, stylusRow.t3];                   % Translation

    % Convert quaternion and translation to 4x4 rigid body transformation matrix
    % The quaternion that comes from the csv file is in the ZYX order
    RStylus = quat2rotm(qStylus);
    TStylus = [RStylus, tStylus(:); 0 0 0 1];

    % Create Stylus struct
    Stylus.name = stylusName;
    Stylus.T = TStylus;

    % Extract rows from 2 to end for Probes data
    probesRows = tempData(2:end, :);
    
    % Initialize an array of structs for Probes
    Probes(height(probesRows), 1) = struct('name', [], 'T', []);
    
    % Loop through each probe row and extract its data
    for j = 1:height(probesRows)

        % Extract each probe
        probesRow = probesRows(j, :);

        % Matlab string can be weird, better check the datatype first
        if iscell(probesRow.name)
            probesName = probesRow.name{1};
        elseif ischar(probesRow.name) || isstring(probesRow.name)
            probesName = probesRow.name;
        else
            error('Unexpected data type for "name" in file: %s', filename);
        end

        % Get the quaternion and the translation
        % The quaternion that comes from the csv file is in the ZYX order
        qProbes = [probesRow.q4, probesRow.q1, probesRow.q2, probesRow.q3];  % Quaternion (q4 is the scalar, i put it in the front using matlab notation)
        tProbes = [probesRow.t1, probesRow.t2, probesRow.t3];                % Translation

        % Convert quaternion and translation to 4x4 rigid body transformation matrix
        RProbes = quat2rotm(qProbes);
        TProbes = [RProbes, tProbes(:); 0 0 0 1];

        % Create Probes struct
        Probes(j).name = probesName;
        Probes(j).T = TProbes;
    end

    % Store timestamp, Stylus, and Probes in the structure array
    allDataStruct(i).Timestamp = timestamp;
    allDataStruct(i).Stylus = Stylus;
    allDataStruct(i).Probes = Probes;

    % Optionally, display the name of the file read
    fprintf('Processed file: %s with timestamp %d\n', filename, timestamp);
end

end