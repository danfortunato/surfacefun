function f = diff(f, n, dim)
%DIFF   Differentiate a SPHEREFUNV.
%   DIFF(F, DIM, N) is the N-th derivative of the SURFACEFUNV F in the
%   dimension DIM:
%      DIM = 1 is the derivative in the x-direction. (default)
%      DIM = 2 is the derivative in the y-direction.
%      DIM = 3 is the derivative in the z-direction.
%
%   DIFF(F, [NX NY NZ]) takes the NX-th derivative of F in the x-direction,
%   the NY-th derivative of F in the y-direction, and the NZ-th derivative
%   of F in the z-direction.
%
%   See also DIFFX, DIFFY, DIFFZ.

% Empty check:
if ( isempty(f) )
    return
end

% Default to first derivative:
if ( nargin < 2 || isempty(n) )
    n = 1;
elseif ( n == 0 )
    return
end

% Default to partial derivative in X:
if ( nargin < 3 )
    dim = 1;
elseif ( numel(dim) ~= 1 )
    error('SURFACEFUNV:diff:dim', 'Dimension should be either 1, 2, or 3.');
end

% Make sure N is in [NX NY NZ] format
if ( isscalar(n) )
    if ( dim == 1 )
        n = [n 0 0];
    elseif ( dim == 2 )
        n = [0 n 0];
    elseif ( dim == 3 )
        n = [0 0 n];
    else
        error('SURFACEFUNV:diff:dim', 'Dimension should be either 1, 2, or 3.');
    end
end

f.components{1} = diff(f.components{1}, n, dim);
f.components{2} = diff(f.components{2}, n, dim);
f.components{3} = diff(f.components{3}, n, dim);

end
