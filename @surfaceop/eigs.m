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
%   EIGS(..., NAME, VALUE) will pass along the name-value arguments to the
%   built-in EIGS routine. See the documentation of EIGS for more details.

if ( ~isInitialized(L) )
    error('SURFACEOP:eigs:notInitialized', ...
        'The surfaceop has not been initialized.');
end

k = 6;
if ( nargin > 1 && isnumeric(varargin{1}) )
    k = varargin{1};
    varargin = varargin(2:end);
end

sigma = 0;
if ( nargin > 2 && isnumeric(varargin{1}) )
    sigma = varargin{1};
    varargin = varargin(2:end);
end

% If a nonzero shift was given then we must rebuild the SURFACEOP:
if ( sigma ~= 0 )
    pdo = L.op;
    pdo.b = pdo.b - sigma;
    rankdef = L.rankdef;
    L = surfaceop(L.domain, pdo);
    L.rankdef = rankdef;
end

if ( ~isBuilt(L) )
    L.build();
end

n = numel(L.domain);
f = surfacefun(L.domain);

    function x = Linv(b)
        f.vec(:) = b;
        if ( L.rankdef )
            f = f - mean(f);
        end
        L = L.updateRHS(f); % Overloaded subsasgn doesn't work here.
        u = L.solve(0);
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
