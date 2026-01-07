function q_conj = quat_conj(q)
%QUAT_CONJ Conjugate.
q_conj = [q(1); -q(2:4)];
end