function q = quat_normalize(q)
%QUAT_NORMALIZE Ensure unit quaternion.
n = norm(q);
if n > 0
    q = q ./ n;
else
    q = [1;0;0;0];
end
end