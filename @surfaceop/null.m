function Z = null(L, tol)
%NULL   Null space of a SURFACEOP.
%   Z = NULL(L) is an orthonormal basis for the null space of the SURFACEOP
%   L. If L has not yet been built (see BUILD()) then NULL(L) will build
%   it. If L has not yet been initialized (see INITIALIZE()) then an error
%   is thrown.
%
%   Z = NULL(L, TOL) uses the tolerance TOL to compute the null space.
%
% See also RANK.

if ( ~isInitialized(L) )
    error('SURFACEOP:null:notInitialized', 'The surfaceop has not been initialized.');
end

if ( ~isBuilt(L) )
    L.build();
end

% Choose a tolerance if none is given:
if ( nargin < 2 )
	tol = 1e-13;
end

% Get the top-level linear system:
A = L.patches{1}.A;

% Compute its SVD and rank:
[~, S, V] = svd(A, 0);
s = diag(S);
r = sum(s/s(1) > tol);
V = V(:,r+1:end);

% If the null space is empty, return []:
if ( size(V, 2) == 0 )
    Z = [];
    return
end

% Propagate each null vector down the tree as a particular solution:
Z = surfacefun;
for k = 1:size(V, 2)
    L.patches{1}.u_part = V(:,k);
    Z(:,k) = L.solve();
end

% Orthonormalize:
Z = orth(Z);

end
