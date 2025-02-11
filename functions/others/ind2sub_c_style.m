function voxelcoordinate = ind2sub_c_style(indices, dimsize)
%IND2SUB_C_STYLE Convert linear indices to subscripts (C-style) for 3D arrays.
%   
%   voxelcoordinate = IND2SUB_C_STYLE(indices, dimsize) returns an N-by-3
%   matrix of coordinates (x, y, z), where:
%       x = mod(indices, dimsize(1));
%       y = mod(floor(indices / dimsize(1)), dimsize(2));
%       z = floor(indices / (dimsize(1) * dimsize(2)));
%
%   The inputs are:
%       indices : A vector of 0-based linear indices, as used in C++.
%       dimsize : A 1x3 vector [nx, ny, nz] specifying the number of 
%                 elements in x, y, and z dimensions, respectively.
%
%   IMPORTANT:
%       - This is *not* the same as MATLAB's built-in ind2sub, 
%         which uses *column-major* (Fortran-style) order.
%       - If you want 1-based indices (standard for MATLAB), 
%         add +1 to x, y, and z.
%
%   Example:
%       % Suppose dimsize = [4, 3, 2]. That is, x=4, y=3, z=2.
%       % Let indices = [0, 1, 5, 6, 23].
%       voxelcoordinate = ind2sub_c_style([0, 1, 5, 6, 23], [4, 3, 2])
%
%       % This should yield (for zero-based):
%       %   [ 0, 0, 0
%       %     1, 0, 0
%       %     1, 1, 0
%       %     2, 1, 0
%       %     3, 2, 1 ]

    % Basic input checks
    if numel(dimsize) ~= 3
        error('dimsize must be a 1x3 vector specifying [nx, ny, nz].');
    end
    if any(dimsize < 1)
        error('All elements of dimsize must be positive integers.');
    end
    
    % Force indices to be a column vector for easy vectorized math
    indices = double(indices(:));    % Ensure double for floor/mod operations
    dimsize = double(dimsize(:));    % Similarly ensure double
    
    % Vectorized computation of C-style subscripts
    x = mod(indices, dimsize(1));
    y = floor(indices ./ dimsize(1));
    z = floor(indices ./ (dimsize(1) * dimsize(2)));
    y = mod(y, dimsize(2));
    
    % Combine into an N-by-3 matrix
    voxelcoordinate = [x, y, z];

    % OPTIONAL: If you need MATLAB-style (1-based) subscripts, uncomment:
    % voxelcoordinate = voxelcoordinate + 1;
end
