function c = merge(a, b, rankdef)
%MERGELR   Merge two patch objects.
%   C = MERGE(A, B) returns a patch C formed by merging the two patches A
%   and B. Typically A and B will be adjacent and any common edges will be
%   eliminated by enforcing continuity and continuity of the derivative
%   across the boundary.

% Parse inputs:
if ( nargin == 0 )
    c = [];
    return
elseif ( nargin == 1 )
    c = a;
    return
elseif ( nargin == 2 )
    rankdef = false;
end

% Compute the indices of intersecting points in a and b.
[i1, i2, s1, s2, flip1, flip2, scl1, scl2, D2N_scl, dom, edges] = intersect(a, b);

% Extract D2N maps:
D2Na = a.D2N; D2Nb = b.D2N;

% Compute new solution operator:
% - The Dirichlet-to-Neumann maps on singular elements have their Jacobians
%   factored out. Therefore, we need to multiply the continuity conditions
%   by the multiplication matrices SCL1 and SCL2. The coordinate maps are
%   such that multiplying D2NA/B by SCL1/2 cancels out, so D2NA/B should
%   only be multiplied by SCL2/1, respectively.
A = -( scl2.*flip1*D2Na(s1,s1)*flip1.' + scl1.*flip2*D2Nb(s2,s2)*flip2.' );
% z = [ scl2.*flip1*D2Na(s1,i1), ...
%       scl1.*flip2*D2Nb(s2,i2), ...
%       scl2.*flip1*D2Na(s1,end) + scl1.*flip2*D2Nb(s2,end) ];
%   |----------------------- rhs -----------------------|
z = [ scl2.*flip1*D2Na(s1,i1) scl1.*flip2*D2Nb(s2,i2) ];
z_part = scl2.*flip1*a.du_part(s1,:) + scl1.*flip2*b.du_part(s2,:);

% Check for a closed surface at the top level
if ( rankdef && isempty(i1) && isempty(i2) )
    % Fix rank deficiency with Leslie's ones matrix trick:
    A = A + ones(size(A));
end

% Store the decomposition for reuse in updateRHS():
dA = decomposition(A);
S = dA \ z;
u_part = dA \ z_part;

if ( isIllConditioned(dA) )
    keyboard
end

% Compute new D2N maps:
Z12 = zeros(numel(i1), numel(i2));
%                                 |--- rhs ----|
% D2N = [ D2Na(i1,i1),  Z12,         D2Na(i1,end) ;
%         Z12.',        D2Nb(i2,i2), D2Nb(i2,end) ] ...
%     + [ D2Na(i1,s1)*flip1.' ; D2Nb(i2,s2)*flip2.' ] * S;
M = [ D2Na(i1,s1)*flip1.' ; D2Nb(i2,s2)*flip2.' ];
D2N = [ D2Na(i1,i1) Z12 ; Z12.' D2Nb(i2,i2) ] + M * S;
du_part = [ a.du_part(i1,:) ; b.du_part(i2,:) ] + M * u_part;

% Construct the new patch:
xyz = [a.xyz(i1,:) ; b.xyz(i2,:)];
id = [a.id ; b.id];
c = surfaceop.parent(dom, id, S, D2N, D2N_scl, u_part, du_part, A, dA, edges, xyz, a, b, ...
    {i1, s1}, {i2, s2}, flip1, flip2, scl1, scl2);

end
