function Q = orth(f, tol)
%ORTH   Array-valued SURFACEFUN orthonormalization.
%   Q = ORTH(F) is an orthonormal basis for the range of the array-valued
%   SURFACEFUN F. That is, Q'*Q = I, the columns of Q span the same space
%   as the columns of F, and the number of columns of Q is the rank of F.
%
%   Q = ORTH(F, TOL) uses the tolerance TOL for orthonormalization.
%
% See also NULL, SVD, RANK, QR.

% Choose a tolerance if none is given:
if ( nargin < 2 )
	tol = 1e-13;
end

% Compute the SVD:
[U, S, ~] = svd(f);
s = diag(S);

% Compute the rank:
r = sum(s/s(1) > tol);

% Select these columns of U:
Q = U(:,1:r);

end
