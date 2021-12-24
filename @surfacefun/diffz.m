function f = diffz(f, n)
%DIFFZ   Differentiate a SURFACEFUN with respect to z.
%   DIFFZ(F) returns a SURFACEFUN representing the derivative of the
%   SURFACEFUN F in the z-direction.
%
%   DIFFZ(F, N) returns an SURFACEFUN representing the N-th derivative of
%   the SURFACEFUN F in the z-direction.
%
%   See also DIFFX, DIFFY, DIFF.

% Default to first derivative:
if ( nargin == 1 )
    n = 1;
end

f = diff(f, n, 3);

end
