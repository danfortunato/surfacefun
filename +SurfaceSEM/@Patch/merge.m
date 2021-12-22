function c = merge(a, b)
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
end

% Compute the indices of intersecting points in a and b.
[i1, i2, s1, s2, dom, edges] = intersect(a, b);

% Extract D2N maps:
D2Na = a.D2N; D2Nb = b.D2N;

% Compute new solution operator:
A = -( D2Na(s1,s1) + D2Nb(s2,s2) );
z = [ D2Na(s1,i1), D2Nb(s2,i2), D2Na(s1,end) + D2Nb(s2,end) ];
%                              |----------- rhs -----------|

% Fix rank deficiency with Leslie's ones matrix trick
if ( rank(A) < size(A,1) )
    A = A + ones(size(A));
end

% Store the decomposition for reuse in updateRHS().
dA = decomposition(A, 'lu');
S = dA \ z;

% Compute new D2N maps:
Z12 = zeros(numel(i1), numel(i2));
%                                 |--- rhs ----|
D2N = [ D2Na(i1,i1),  Z12,         D2Na(i1,end) ;
        Z12.',        D2Nb(i2,i2), D2Nb(i2,end) ] ...
    + [ D2Na(i1,s1) ; D2Nb(i2,s2) ] * S;

% Construct the new patch:
xyz = [a.xyz(i1,:) ; b.xyz(i2,:)];
c = SurfaceSEM.Parent(dom, S, D2N, dA, edges, xyz, a, b, {i1, s1}, {i2, s2});

end
