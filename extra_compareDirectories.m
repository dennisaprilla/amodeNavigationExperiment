%% Define the paths to Directory A and Directory B

root_dir = "D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_a\data\trial_0006_Session2_02\snapshot";
dirA = fullfile(root_dir, "A_T_MID");  % Replace with the actual path for Directory A
dirB = fullfile(root_dir, "A_T_MID_reobserved");  % Replace with the actual path for Directory B

%% Get file lists for both directories (ignoring subdirectories)
filesA = dir(fullfile(dirA, '*'));
filesA = filesA(~[filesA.isdir]);  % Remove directories

filesB = dir(fullfile(dirB, '*'));
filesB = filesB(~[filesB.isdir]);  % Remove directories

%% Extract unique identifier numbers from Directory A
IDs_A = {};  % Initialize as a cell array of strings
for i = 1:length(filesA)
    % Use regexp to extract digits at the beginning of the filename.
    % This assumes the identifier is a sequence of digits at the start.
    tokens = regexp(filesA(i).name, '^(\d+)', 'tokens');
    if ~isempty(tokens)
        IDs_A{end+1} = tokens{1}{1};
    end
end
IDs_A = unique(IDs_A);  % Ensure we have unique identifiers

%% Extract unique identifier numbers from Directory B
IDs_B = {};
for i = 1:length(filesB)
    tokens = regexp(filesB(i).name, '^(\d+)', 'tokens');
    if ~isempty(tokens)
        IDs_B{end+1} = tokens{1}{1};
    end
end
IDs_B = unique(IDs_B);

%% Compare the identifier lists: find IDs in A that are missing in B
missingIDs = setdiff(IDs_A, IDs_B);

%% Display the result
if isempty(missingIDs)
    fprintf('No sets were deleted. All identifier sets from Directory A are present in Directory B.\n');
else
    fprintf('The following identifier numbers are missing in Directory B:\n');
    for i = 1:length(missingIDs)
        fprintf('%s\n', missingIDs{i});
    end
end
