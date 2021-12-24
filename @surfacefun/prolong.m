function out = prolong(f, n)
%PROLONG   Prolong a SURFACEFUN.
%   PROLONG(F, N) returns a SURFACEFUN representing the same function as F,
%   but using an N x N discretization on each element.

if ( nargin == 1 )
    out = f;
    return
end

out = f;
for k = 1:length(f)
    [nx, ny] = size(f.vals{k});
    ny = min(ny, n);
    nx = min(nx, n);
    coeffs = chebtech2.vals2coeffs(chebtech2.vals2coeffs(f.vals{k}).').';
    U = zeros(n);
    U(1:ny,1:nx) = coeffs(1:ny,1:nx);
    out.vals{k} = chebtech2.coeffs2vals(chebtech2.coeffs2vals(U).').';
end

end
