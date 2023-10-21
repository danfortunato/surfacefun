function c = merge_ItI(a, b, rankdef)
%MERGE   Merge two patch objects.
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
[i1, i2, s1, s2, flip1, flip2, scl1, scl2, ItI_scl, dom, edges] = intersect(a, b);

% Extract ItI maps:
ItIa = a.BtB;
ItIb = b.BtB;

% Compute new solution operator:
R11a = ItIa(i1,i1);               R22b = ItIb(i2,i2);
R13a = ItIa(i1,s1)*flip1.';       R23b = ItIb(i2,s2)*flip2.';
R31a = flip1*ItIa(s1,i1);         R32b = flip2*ItIb(s2,i2);
R33a = flip1*ItIa(s1,s1)*flip1.'; R33b = flip2*ItIb(s2,s2)*flip2.';
I = eye(size(R33a));

W = I - R33b*R33a;

% Check for a closed surface at the top level:
if ( rankdef && isempty(i1) && isempty(i2) )
    w = a.w(s1);
    W = W + w*w';
end

W = inv(W); %#ok<*MINV>

h1a = a.du_part(i1,:);
h2b = b.du_part(i2,:);
h3a = flip1*a.du_part(s1,:);
h3b = flip2*b.du_part(s2,:);

W_R33b_R31a = W*R33b*R31a;
W_R32b = W*R32b;

% Compute new solution operators:
Sa = [ W_R33b_R31a, -W_R32b ];
Sb = [ -R31a - R33a*W_R33b_R31a, R33a*W_R32b ];

u_part_a = W*(R33b*h3a - h3b);
u_part_b = -(I + R33a*W*R33b)*h3a + R33a*(W*h3b);

% Compute new D2N map:
ItI = [  R11a + R13a*W_R33b_R31a        -R13a*W_R32b             ;
        -R23b*(R31a + R33a*W_R33b_R31a)  R22b + R23b*R33a*W_R32b ];

du_part = [ h1a + R13a*u_part_a ;
            h2b + R23b*u_part_b ];

% Store the decomposition for reuse in updateRHS():
A = I - R33b*R33a;
dA = decomposition(A);

% Construct the new patch:
xyz = [a.xyz(i1,:) ; b.xyz(i2,:)];
w = [a.w(i1) ; b.w(i2)];
id = [a.id ; b.id];
c = surfaceop.parent(dom, id, {Sa, Sb}, ItI, ItI_scl, {u_part_a, u_part_b}, du_part, A, dA, ...
    edges, xyz, w, a, b, {i1, s1}, {i2, s2}, flip1, flip2, scl1, scl2);

end
