function P = updateRHS(P, rhs)
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
a = updateRHS(P.child1, rhs(1:n1));
b = updateRHS(P.child2, rhs(n1+1:n1+n2));

i1 = P.idx1{1};
s1 = P.idx1{2};
i2 = P.idx2{1};
s2 = P.idx2{2};
flip1 = P.flip1;
flip2 = P.flip2;
scl1 = P.scl1;
scl2 = P.scl2;

% Extract D2N maps:
D2Na = a.D2N; D2Nb = b.D2N;
% and discard from children
% a.D2N = []; b.D2N = [];

% Compute new solution operator:
%S = P.dA \ (scl2.*(flip1*D2Na(s1,end)) + scl1.*(flip2*D2Nb(s2,end)));
S = P.dA \ (scl2.*(flip1*a.du_part(s1,:)) + scl1.*(flip2*b.du_part(s2,:)));
%          |------------------------- rhs -------------------------|

% TODO: Also scale this???
% Compute new D2N maps:
%      |--- rhs ----|
% D2N = [ D2Na(i1,end) ;
%         D2Nb(i2,end) ] ...
%     + [ D2Na(i1,s1)*flip1.' ; D2Nb(i2,s2)*flip2.' ] * S;
%      |--- rhs ----|
D2N = [ a.du_part(i1,:) ;
        b.du_part(i2,:) ] ...
    + [ D2Na(i1,s1)*flip1.' ; D2Nb(i2,s2)*flip2.' ] * S;

%P.S(:,end) = S;
%P.D2N(:,end) = D2N;
P.u_part = S;
P.du_part = D2N;

P.child1 = a;
P.child2 = b;

end
