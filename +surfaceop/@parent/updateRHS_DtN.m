function P = updateRHS_DtN(P, rhs)
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
a = updateRHS_DtN(P.child1, rhs(1:n1,:));
b = updateRHS_DtN(P.child2, rhs(n1+1:n1+n2,:));

i1 = P.idx1{1};
s1 = P.idx1{2};
i2 = P.idx2{1};
s2 = P.idx2{2};
flip1 = P.flip1;
flip2 = P.flip2;
scl1 = P.scl1;
scl2 = P.scl2;

% Extract D2N maps:
DtNa = a.BtB; DtNb = b.BtB;
% and discard from children
% a.D2N = []; b.D2N = [];

% Compute new solution operator:
P.u_part = P.dA \ (scl2.*(flip1*a.du_part(s1,:)) + scl1.*(flip2*b.du_part(s2,:)));

% Compute new D2N map:
% TODO: Also scale this???
P.du_part = [ a.du_part(i1,:) ; b.du_part(i2,:) ] + ...
            [ DtNa(i1,s1)*flip1.' ; DtNb(i2,s2)*flip2.' ] * P.u_part;

P.child1 = a;
P.child2 = b;

end
