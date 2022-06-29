function r = rank(f, tol)
%RANK   Rank of an array-valued SURFACEFUN.
%   RANK(F) produces an estimate of the number of linearly independent
%   columns of the array-valued SURFACEFUN F.
%
%   RANK(F, TOL) is the number of relative singular values of F greater
%   than TOL.
%
% See also SVD, QR, ORTH.

% Choose a tolerance if none is given:
if ( nargin < 2 )
	tol = 1e-13;
end

% Compute the singular values:
s = svd(f);

% Compute the rank:
r = sum(s/s(1) > tol);

end
