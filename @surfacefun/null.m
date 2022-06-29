function Z = null(f, tol)
%NULL   Null space of an array-valued SURFACEFUN.
%   Z = NULL(F) is an orthonormal basis for the null space of the
%   array-valued SURFACEFUN F.
%
%   Z = NULL(F, TOL) uses the tolerance TOL to compute the null space.
%
% See also ORTH, SVD, RANK, QR.

% Choose a tolerance if none is given:
if ( nargin < 2 )
	tol = 1e-13;
end

% Compute the SVD:
[~, S, V] = svd(f);
s = diag(S);

% Compute the rank:
r = sum(s/s(1) > tol);

% The null vectors are the final columns of V:
Z = V(:,r+1:end);

end
