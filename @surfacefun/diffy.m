function f = diffy(f, n)
%DIFFY   Differentiate a SURFACEFUN with respect to y.
%   DIFFY(F) returns a SURFACEFUN representing the derivative of the
%   SURFACEFUN F in the y-direction.
%
%   DIFFY(F, N) returns a SURFACEFUN representing the N-th derivative of
%   the SURFACEFUN F in the y-direction.
%
%   See also DIFFX, DIFFZ, DIFF.

% Default to first derivative:
if ( nargin == 1 )
    n = 1;
end

f = diff(f, n, 2);

end
