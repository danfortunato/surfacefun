function f = diff(f, n, dim)
%DIFF   Differentiate a SURFACEFUN.
%   DIFF(F, N, DIM) is the N-th derivative of the SURFACEFUN F in the
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
    error('SURFACEFUN:diff:dim', 'Dimension should be either 1, 2, or 3.');
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
        error('SURFACEFUN:diff:dim', 'Dimension should be either 1, 2, or 3.');
    end
end
nx = n(1);
ny = n(2);
nz = n(3);

for k = 1:length(f)
    vals = f.vals{k};
    coeffs = chebtech2.vals2coeffs( chebtech2.vals2coeffs(vals).' ).';
    for m = 1:nx, coeffs = mappedDiff(coeffs, f.domain, k, 1); end
    for m = 1:ny, coeffs = mappedDiff(coeffs, f.domain, k, 2); end
    for m = 1:nz, coeffs = mappedDiff(coeffs, f.domain, k, 3); end
    vals = chebtech2.coeffs2vals( chebtech2.coeffs2vals(coeffs).' ).';
    f.vals{k} = vals;
end

end

function f = mappedDiff(f, dom, k, dim)
%MAPPEDDIFF   Mapped derivative of a matrix of Chebyshev coefficients.
%   MAPPEDDIFF(F, DOM, K, DIM) returns the matrix of Chebyshev coefficients
%   representing the mapped derivative of F in the dimension DIM. The
%   mapping is defined by the K-th coordinate transformations in DOM.

    % Compute derivatives on [-1,1]^2
    [nv, nu] = size(f);
    dfdu = [ cdiff(f.').', zeros(nv, 1) ];
    dfdv = [ cdiff(f);     zeros(1, nu) ];
    
    %D = diffmat(n);    
    %dfdu2 = f * D.';
    %dfdv2 = D * f;

    % Convert to values to multiply by Jacobian factors
    dfdu = chebtech2.coeffs2vals( chebtech2.coeffs2vals(dfdu).' ).';
    dfdv = chebtech2.coeffs2vals( chebtech2.coeffs2vals(dfdv).' ).';

    % Get Jacobian factors for the specified dimension
    if ( dim == 1 )
        du = dom.ux{k};
        dv = dom.vx{k};
    elseif ( dim == 2 )
        du = dom.uy{k};
        dv = dom.vy{k};
    elseif ( dim == 3 )
        du = dom.uz{k};
        dv = dom.vz{k};
    else
        error('SURFACEFUN:mappedDiff:dim', 'Dimension should be either 1 or 2.');
    end

    % Compute df/dx = (du/dx) df/du + (dv/dx) df/dv
    vals = du .* dfdu + dv .* dfdv;
    if ( dom.singular(k) )
        % Computing derivatives will be incorrect on these patches, since
        % we do not divide by the Jacobian. We could perhaps do something
        % like:
        %
        %    idx = abs(dom.J{k}) > 1e-4;
        %    vals(idx(:)) = vals(idx) ./ dom.J{k}(idx);
        %
        % But for now, we live with the incorrect answer on these patches.
    end

    %Dx = ux(:).*Du + vx(:).*Dv;
    %Dy = uy(:).*Du + vy(:).*Dv;
    %Dz = uz(:).*Du + vz(:).*Dv;

    % Convert back to coefficients
    f = chebtech2.vals2coeffs( chebtech2.vals2coeffs(vals).' ).';

end

function dC = cdiff(C)
%CDIFF   Recurrence relation for coefficients of derivative.
%   CDIFF(C) returns the matrix of Chebyshev coefficients whose columns are
%   the derivatives of the columns of C.

[n, m] = size(C);
dC = zeros(n-1, m);                        % Initialize vector {c_r}
w = repmat(2*(1:n-1)', 1, m);
v = w.*C(2:end,:);                         % Temporal vector
dC(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:), 1); % Compute c_{n-2}, c_{n-4}, ...
dC(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:), 1); % Compute c_{n-3}, c_{n-5}, ...
dC(1,:) = .5*dC(1,:);                      % Adjust the value for c_0

end
