function thresholdSliderExample(volumeData, headerInfo, ax)
    % thresholdSliderExample(volumeData, headerInfo, ax)
    %
    % Displays a slider that controls a threshold, and a scatter plot of
    % the voxel points in 'volumeData' that exceed that threshold.  All of
    % this is done within an existing figure/axes specified by 'ax'.
    %
    % Inputs:
    %    volumeData  : 1D or 3D array containing your volume data
    %    headerInfo  : struct with fields DimSize, Offset, ElementSpacing, 
    %                  and TransformMatrix
    %    ax          : handle to the axes where the scatter plot should appear
    %
    % Requires:
    %    The function ind2sub_c_style (or replace with MATLAB's ind2sub).
    %
    % ------------
    % Example:
    %   fig = figure;
    %   ax  = axes('Parent',fig);
    %   volumeData = rand(50,50,50) * 200; % dummy data
    %   headerInfo.DimSize = [50,50,50];
    %   headerInfo.Offset  = [0,0,0];
    %   headerInfo.ElementSpacing  = [1,1,1];
    %   headerInfo.TransformMatrix = reshape(eye(3),1,[]);
    %   thresholdSliderExample(volumeData, headerInfo, ax);

    % Validate inputs
    if nargin < 3
        error('Need to provide volumeData, headerInfo, and axes handle as inputs.');
    end
    if ~ishandle(ax) || ~strcmpi(get(ax, 'type'), 'axes')
        error('The third input must be a valid axes handle.');
    end
    
    % Get the figure that contains these axes
    fig = ancestor(ax, 'figure');

    % Store the data in the figure’s UserData (so the slider callback can see it)
    fig.UserData.volumeData  = volumeData;
    fig.UserData.headerInfo  = headerInfo;

    % Initial threshold
    initialThreshold = 100;

    % Create a slider in the figure (parent is the figure, not the axes)
    % Adjust 'Position' and min/max values as appropriate for your data
    thresholdSlider = uicontrol('Parent', fig, ...
        'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.10 0.90 0.80 0.05], ...  % [left bottom width height]
        'Min', 0, ...
        'Max', max(volumeData(:)), ...
        'Value', initialThreshold, ...
        'Callback', @(src,~) updateScatter(src, ax, fig));

    % Initial plot
    updateScatter(thresholdSlider, ax, fig);
end


function updateScatter(src, ax, fig)
    % Update the scatter plot when the slider changes

    % Retrieve data from the figure
    volumeData = fig.UserData.volumeData;
    headerInfo = fig.UserData.headerInfo;

    % Current slider threshold
    threshold = src.Value;

    % Find indices above threshold
    idx_overthresh = find(volumeData > threshold);

    % Convert linear indices to voxel coordinates
    % (Replace 'ind2sub_c_style' with your actual function or 'ind2sub')
    voxelcoordinate = ind2sub_c_style(idx_overthresh, headerInfo.DimSize);

    % Convert to homogeneous coordinates
    voxelcoordinate_homogeneous = [voxelcoordinate'; ones(1, size(voxelcoordinate,1))];

    % Build transformation matrix from header info
    init_t = headerInfo.Offset;                           % 1x3
    init_s = diag(headerInfo.ElementSpacing);             % 3x3
    init_R = reshape(headerInfo.TransformMatrix, [3, 3]); % 3x3
    init_A = [init_R .* init_s, init_t'; 0 0 0 1];        % 4x4

    % Apply the transformation
    voxelcoordinate_homogeneous = init_A * voxelcoordinate_homogeneous;

    % Check if scatter already exists
    hScatter = findobj(ax, 'Type', 'Scatter', 'Tag', 'bone_surface');
    if isempty(hScatter)
        % Create a new scatter if it does not exist
        scatter3(ax, ...
            voxelcoordinate_homogeneous(1,:), ...
            voxelcoordinate_homogeneous(2,:), ...
            voxelcoordinate_homogeneous(3,:), ...
            10, 'filled', ...
            'Tag', 'bone_surface');
    else
        % Update existing scatter
        set(hScatter, ...
            'XData', voxelcoordinate_homogeneous(1,:), ...
            'YData', voxelcoordinate_homogeneous(2,:), ...
            'ZData', voxelcoordinate_homogeneous(3,:));
    end

    % Force MATLAB to redraw
    drawnow;
end
