function f = diffx(f, n)
%DIFFX   Differentiate a SURFACEFUN with respect to x.
%   DIFFX(F) returns a SURFACEFUN representing the derivative of the
%   SURFACEFUN F in the x-direction.
%
%   DIFFX(F, N) returns an SURFACEFUN representing the N-th derivative of
%   the SURFACEFUN F in the x-direction.
%
%   See also DIFFY, DIFFZ, DIFF.

% Default to first derivative:
if ( nargin == 1 )
    n = 1;
end

f = diff(f, n, 1);

end
