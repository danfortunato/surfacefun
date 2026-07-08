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

dom = f(1).domain;
if ( ~all(dom.ptype == dom.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

switch ( dom.ptype(1) )
    case 'tri'
        [K, Ku, Kv] = koornwinder(order(dom));
        Du = Ku / K;
        Dv = Kv / K;
        du = @(x) Du * x;
        dv = @(x) Dv * x;
    case 'quad'
        %[nv, nu] = size(f(1).vals{1});
        [nu, nv] = size(dom);
        Du = diffmat(nu);
        Dv = diffmat(nv);
        du = @(x) x * Du.';
        dv = @(x) Dv * x;
end

% TODO: This can be vectorized across multiple functions.
nf = builtin('numel', f);
for j = 1:nf
    for k = 1:length(f(j))
        vals = f(j).vals{k};
        for m = 1:nx, vals = mappedVDiff(vals, f(j).domain, k, 1, du, dv); end
        for m = 1:ny, vals = mappedVDiff(vals, f(j).domain, k, 2, du, dv); end
        for m = 1:nz, vals = mappedVDiff(vals, f(j).domain, k, 3, du, dv); end
        f(j).vals{k} = vals;
    end
end

end

function f = mappedVDiff(f, dom, k, dim, du, dv)

    dfdu = du(f);
    dfdv = dv(f);
    %dfdu = f * Du.';
    %dfdv = Dv * f;

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
        error('SURFACEFUN:mappedVDiff:dim', 'Dimension should be either 1 or 2.');
    end

    % Compute df/dx = (du/dx) df/du + (dv/dx) df/dv
    f = du .* dfdu + dv .* dfdv;
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

end
