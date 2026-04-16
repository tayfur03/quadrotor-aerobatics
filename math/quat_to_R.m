function R = quat_to_R(q)
%QUAT_TO_R Convert quaternion to rotation matrix (body->inertial).
% q = [qw; qx; qy; qz]

qw = q(1); qx = q(2); qy = q(3); qz = q(4);

R = [...
    1-2*(qy^2+qz^2),  2*(qx*qy - qw*qz),  2*(qx*qz + qw*qy);
    2*(qx*qy + qw*qz), 1-2*(qx^2+qz^2),  2*(qy*qz - qw*qx);
    2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1-2*(qx^2+qy^2)];
end