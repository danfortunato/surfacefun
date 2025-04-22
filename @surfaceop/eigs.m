function varargout = eigs(L, varargin)
%EIGS   Find a few eigenvalues and eigenfunctions of a SURFACEOP.
%   D = EIGS(L) returns a vector of 6 eigenvalues of the SURFACEOP L. EIGS
%   will attempt to return the eigenvalues corresponding to the smoothest
%   eigenfunctions.
%
%   [V, D] = EIGS(L) returns a diagonal matrix D of L's 6 smoothest
%   eigenvalues and an array-valued SURFACEFUN V of the corresponding
%   eigenfunctions.
%
%   [V, D, FLAG] = EIGS(L) also returns a convergence flag. If FLAG is 0
%   then all the eigenvalues converged; otherwise not all converged.
%
%   EIGS(L, K) returns the K smoothest eigenvalues.
%
%   EIGS(L, K, SIGMA) returns the K smoothest eigenvalues closest to the
%   scalar SIGMA, which may be real or complex, including 0.
%
%   EIGS(..., 'dirichlet') or EIGS(..., 'neumann') computes Dirichlet or
%   Neumann eigenvalues, respectively.
%
%   EIGS(..., NAME, VALUE) will pass along the name-value arguments to the
%   built-in EIGS routine. See the documentation of EIGS for more details.

if ( ~isInitialized(L) )
    error('SURFACEOP:eigs:notInitialized', ...
        'The surfaceop has not been initialized.');
end

k = 6;
sigma = 0;
bc = 'dirichlet';
parsed = [];
parsed.k     = false;
parsed.sigma = false;
parsed.bc    = false;

if ( nargin > 1 )
    if ( isnumeric(varargin{1}) )
        k = varargin{1};
        varargin = varargin(2:end);
        parsed.k = true;
    elseif ( (isstring(varargin{1}) || ischar(varargin{1})) && ...
             any(strcmpi(varargin{1}, {'dirichlet', 'neumann'})) )
        bc = lower(varargin{1});
        varargin = varargin(2:end);
        parsed.bc = true;
    end
end

if ( nargin > 2 )
    if ( parsed.bc && ~parsed.k && isnumeric(varargin{1}) )
        k = varargin{1};
        varargin = varargin(2:end);
        parsed.k = true;
    elseif ( parsed.k && isnumeric(varargin{1}) )
        sigma = varargin{1};
        varargin = varargin(2:end);
        parsed.sigma = true;
    elseif ( parsed.k && (isstring(varargin{1}) || ischar(varargin{1})) && ...
             any(strcmpi(varargin{1}, {'dirichlet', 'neumann'})) )
        bc = lower(varargin{1});
        varargin = varargin(2:end);
        parsed.bc = true;
    end
end

if ( nargin > 3 )
    if ( parsed.k && parsed.bc && ~parsed.sigma && isnumeric(varargin{1}) )
        sigma = varargin{1};
        varargin = varargin(2:end);
        parsed.sigma = true;
    elseif ( parsed.k && parsed.sigma && ~parsed.bc && ...
            (isstring(varargin{1}) || ischar(varargin{1})) && ...
            any(strcmpi(varargin{1}, {'dirichlet', 'neumann'})) )
        bc = lower(varargin{1});
        varargin = varargin(2:end);
        parsed.bc = true;
    end
end

% If a nonzero shift was given then we must rebuild the SURFACEOP:
if ( sigma ~= 0 )
    pdo = L.op;
    pdo.b = pdo.b - sigma;
    L = surfaceop(L.domain, pdo);
end

if ( ~isBuilt(L) )
    L.build();
end

if ( strcmpi(bc, 'neumann') )
    % Map homogeneous Neumann data to equivalent Dirichlet data
    DtN = L.DtN();
    nb = size(DtN, 1);
    W = ones(nb);
    dA = decomposition(DtN + W, 'CheckCondition', false);
    dir_bc = @(x) dA \ -x;
else
    dir_bc = @(x) 0;
end

n = numel(L.domain);
f = surfacefun(L.domain);

    function x = Linv(b)
        f.vec(:) = b;
        if ( L.rankdef )
            f = f - mean(f);
        end
        L = L.updateRHS(f); % Overloaded subsasgn doesn't work here.
        g = dir_bc(L.patches{1}.du_part);
        u = L.solve(g);
        x = u.vec();
    end

% Call MATLAB's built-in EIGS routine:
[VV, D, flag] = eigs(@Linv, n, k, sigma, varargin{:});

% Convert the eigenvectors to an array-valued SURFACEFUN:
V = repmat(surfacefun(L.domain), 1, k);
V.vec(:) = VV;

if ( L.rankdef )
    % Make sure the eigenvectors are mean-zero
    V = V - mean(V);
end

if ( nargout < 2 )
    varargout = {diag(D)};
else
    varargout = {V, D, flag};
end

end
