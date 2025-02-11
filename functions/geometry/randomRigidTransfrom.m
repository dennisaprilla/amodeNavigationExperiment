function T = randomRigidTransfrom(angleLimit, translationLimit)
    % Generate random small rotations within the range [-10, 10] degrees
    rx = deg2rad(angleLimit * (2 * rand - 1)); % Rotation about x-axis
    ry = deg2rad(angleLimit * (2 * rand - 1)); % Rotation about y-axis
    rz = deg2rad(angleLimit * (2 * rand - 1)); % Rotation about z-axis

    % Construct individual rotation matrices
    Rx = [1  0       0;
          0  cos(rx) -sin(rx);
          0  sin(rx)  cos(rx)];

    Ry = [cos(ry)  0  sin(ry);
          0        1  0;
         -sin(ry)  0  cos(ry)];

    Rz = [cos(rz) -sin(rz)  0;
          sin(rz)  cos(rz)  0;
          0        0        1];

    % Combine into a single rotation matrix
    R = Rz * Ry * Rx; % Order matters: Apply X, then Y, then Z rotation

    % Generate a random translation vector within [-10, 10]
    t = translationLimit * (2 * rand(3,1) - 1);

    % Construct the 4x4 transformation matrix
    T = [R, t; 0 0 0 1];
end
