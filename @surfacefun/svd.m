function [U, S, V] = svd(f)
%SVD   Singular value decomposition of an array-valued SURFACEFUN.
%   [U, S, V] = SVD(F), where F is an array-valued SURFACEFUN with N
%   columns, produces an N x N diagonal matrix S with nonnegative diagonal
%   elements in nonincreasing order, an array-valued SURFACEFUN U with N
%   orthonormal columns, and an N x N unitary matrix V such that F = U*S*V'.
%
%   S = SVD(F) returns a vector containing the singular values of F.

% Call surfacefun QR:
[Q, R] = qr(f);

% Call discrete SVD:
[U, S, V] = svd(R, 0);

% Make U a surfacefun:
U = Q*U;

% Output only the singular values:
if ( nargout < 2 )
    U = diag(S);
end

end
