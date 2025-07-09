function [p] = MNQuatMultiply(q, r)
%QUATMULTIPLY Multiply two quaternions.
%   P = QUATMULTIPLY(Q, R) multiplies quaternions Q and R, returning their
%   product P.
    if size(q, 1) ~= 4 || size(r, 1) ~= 4
      error('Expecting quaternions as rows.');
    elseif size(q, 2) ~= size(r, 2)
      error('Number of quaternions don''t match.');
    end
    q0 = q(1);
    qv = q(2:4);
    r0 = r(1);
    rv = r(2:4);
    
    p = [q0*r0 - dot(qv, rv); q0*rv + r0*qv + cross(qv, rv)];
    p = p / norm(p);
end