function P = updateRHS_ItI(P, rhs)
%UPDATERHS   Update RHS of a SURFACEOP.PARENT object.
%   P = UPDATERHS(P, RHS) replaces the existing RHS of an initialized
%   SURFACEOP.PARENT object P with that given in RHS, which must be a cell
%   array containing tensor-product Chebyshev values for each patch.

if ( ~iscell(rhs) )
    error('SURFACEOP:PARENT:updateRHS:format', 'RHS must be a cell array.');
end

% We are updating the RHS from a cell array of values.
% Get the number of patches in each child:
n1 = P.child1.len;
n2 = P.child2.len;
% Update RHS of children:
a = updateRHS_ItI(P.child1, rhs(1:n1,:));
b = updateRHS_ItI(P.child2, rhs(n1+1:n1+n2,:));

i1 = P.idx1{1};
s1 = P.idx1{2};
i2 = P.idx2{1};
s2 = P.idx2{2};
flip1 = P.flip1;
flip2 = P.flip2;
scl1 = P.scl1;
scl2 = P.scl2;

% Extract D2N maps:
ItIa = a.BtB; ItIb = b.BtB;
% and discard from children
% a.BtB = []; b.BtB = [];

h1a = a.du_part(i1,:);
h2b = b.du_part(i2,:);
h3a = flip1*a.du_part(s1,:);
h3b = flip2*b.du_part(s2,:);

R13a = ItIa(i1,s1)*flip1.';       R23b = ItIb(i2,s2)*flip2.';
R33a = flip1*ItIa(s1,s1)*flip1.'; R33b = flip2*ItIb(s2,s2)*flip2.';
I = eye(size(R33a));

u_part_a = P.dA \ (R33b*h3a - h3b);
u_part_b = -(I + R33a*(P.dA\R33b))*h3a + R33a*(P.dA\h3b);
P.u_part = {u_part_a, u_part_b};
P.du_part = [ h1a + R13a*u_part_a ;
              h2b + R23b*u_part_b ];

% Compute new solution operator:
%P.u_part = P.dA \ (scl2.*(flip1*a.du_part(s1,:)) + scl1.*(flip2*b.du_part(s2,:)));

% Compute new D2N map:
% TODO: Also scale this???
%P.du_part = [ a.du_part(i1,:) ; b.du_part(i2,:) ] + ...
%            [ ItIa(i1,s1)*flip1.' ; ItIb(i2,s2)*flip2.' ] * P.u_part;

P.child1 = a;
P.child2 = b;

end
