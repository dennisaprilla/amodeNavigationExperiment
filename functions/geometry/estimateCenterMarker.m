% Function to fit a sphere and calculate its center
function center = estimateCenterMarker(points)
    % Formulate the matrix for least squares fitting
    A = [2 * points, ones(size(points, 1), 1)];
    b = sum(points.^2, 2);
    % Solve the linear system
    x = A \ b;
    % Extract sphere center
    center = x(1:3);
end
