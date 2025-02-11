function rigidbodyT = estimateRBfrom3Points(points)
    % Estimate the rigid body transformation (rotation + translation)
    % from three 3D points.
    %
    % INPUT:
    %   points - A 3x3 matrix where each column represents a 3D point
    %
    % OUTPUT:
    %   rigidbodyT - A 4x4 homogeneous transformation matrix
    
    % Validate input
    if size(points, 1) ~= 3 || size(points, 2) ~= 3
        error('Input must be a 3x3 matrix, with each column representing a 3D point.');
    end
    
    % Extract the first three points
    p1 = points(:, 1);
    p2 = points(:, 2);
    p3 = points(:, 3);

    % Step 1: Use the first marker as the new origin
    centroid = p1;

    % Step 2: Translate the points so that the centroid is at the origin
    p1_prime = p1 - centroid;
    p2_prime = p2 - centroid;
    p3_prime = p3 - centroid;

    % Step 3: Define the first vector (V1) as p2' - p1'
    v1 = p2_prime - p1_prime;
    v1 = v1 / norm(v1); % Normalize V1

    % Step 4: Define the second vector (V2) as p3' - p1'
    v2 = p3_prime - p1_prime;

    % Orthogonalize V2 by subtracting the projection of V2 onto V1
    v2_proj = (dot(v1, v2) / dot(v1, v1)) * v1;
    v2_orthogonal = v2 - v2_proj;
    v2_orthogonal = v2_orthogonal / norm(v2_orthogonal); % Normalize the orthogonalized V2

    % Step 5: Compute the third orthogonal vector (V3) using the cross product
    v3 = cross(v1, v2_orthogonal);
    v3 = v3 / norm(v3); % Normalize V3

    % Step 6: Construct the rotation matrix from the orthonormal vectors
    rotation = [v1, v2_orthogonal, v3];

    % Step 7a: Apply SVD to ensure the matrix is orthogonal
    [U, ~, V] = svd(rotation);
    rotation = U * V';

    % Step 8: Create the transformation matrix (4x4 homogeneous matrix)
    rigidbodyT = eye(4); % Initialize as identity matrix
    rigidbodyT(1:3, 1:3) = rotation; % Set the rotation part
    rigidbodyT(1:3, 4) = centroid;   % Set the translation part
end