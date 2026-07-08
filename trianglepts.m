function [x, y, w, t2v] = trianglepts(n, opts)
%TRIANGLEPTS   Generate nodes on the triangle.
%   [X, Y] = TRIANGLEPTS(N) generates an order-N set of nodes (X,Y) in the
%   (0,0)-(1,0)-(0,1) unit triangle using Isaac's recursive algorithm [1]. X and
%   Y are vectors of length N*(N+1)/2. By default, the starting 1D node set is
%   taken to be the second-kind Chebyshev nodes (i.e., Chebyshev nodes with
%   endpoints). The resulting nodes (X,Y) are fully symmetric in the triangle
%   and match the starting 1D node set on the edges of the triangle.
%
%   [X, Y] = TRIANGLEPTS(..., type='cheb2') is equivalent to the above. Other
%   possible values for the node type include:
%
%       First-kind Chebyshev:  'cheb1' 'gc'  'chebyshev1'
%       Second-kind Chebyshev: 'cheb2' 'lgc' 'chebyshev2' 'cheb'
%       Gauss-Legendre:        'leg'   'gl'  'legendre'
%       Gauss-Lobatto:         'lob'   'lgl' 'lobatto'
%       Equally spaced nodes:  'equi'  'uni' 'lin' 'equispaced' 'uniform' 'linspace'
%
%   [X, Y, W] = TRIANGLEPTS(N) also returns a vector of order-N interpolatory
%   quadrature weights W for integration on the unit triangle.
%
%   [X, Y, W, T2V] = TRIANGLEPTS(N) also returns the (N-1)^2 x 3 matrix T2V
%   encoding triangle-to-vertex connectivity of the Delaunay triangulation of
%   the nodes (X,Y). The vertices of the i-th triangle are the elements of X and
%   Y with indices given by T2V(i,:).
%
%   [X, Y, ~, T2V] = TRIANGLEPTS(..., weights=false) skips computation of the
%   quadrature weights W.
%
%   References:
%
%      [1] Tobin Isaac, "Recursive, parameter-free, explicitly defined
%          interpolation nodes for simplices", SIAM J. Sci. Comput., 42 (2020),
%          pp. A4046-A4062.

arguments
    n
    opts.dim    = 2
    opts.type   = 'cheb2'
    opts.domain = 'unit'
    opts.weights {mustBeNumericOrLogical} = true
end

switch lower(opts.type)
    case {'leg', 'legendre', 'gl'}
        family = @(n) legpts(n, [0 1]);
    case {'lob', 'lobatto', 'lgl'}
        family = @(n) lobpts_(n, [0 1]);
    case {'cheb1', 'gc', 'chebyshev1'}
        family = @(n) chebpts(n, [0 1], 1);
    case {'cheb2', 'lgc', 'chebyshev2', 'cheb'}
        family = @(n) chebpts(n, [0 1], 2);
    case {'equi', 'equispaced', 'uni', 'uniform', 'lin', 'linspace'}
        family = @(n) linspace(0, 1, n).';
    otherwise
        error('Unknown node set.');
end

dim = opts.dim;

% The number of polynomials up to degree n in dim dimensions
N = nchoosek(n+dim-1, dim);
x = zeros(N, dim+1);
i = 1;
for idx = tuples(dim+1, n-1)
    x(i,:) = recursive(dim+1, n-1, idx+1, family);
    i = i+1;
end
x = from_unit(to_unit(x, 'barycentric'), opts.domain);

y = x(:,2);
x = x(:,1);

w = [];
if ( nargout > 2 && opts.weights )
    % Compute interpolatory quadrature weights
    K = koornwinder(n-1, x, y);
    g = zeros(1, size(K,2));
    g(1) = sqrt(2)/2;
    w = g / K;
    w = w(:);
end

if ( nargout > 3 )
    % Compute triangle-to-vertex connectivity
    t2v = trilattice(n);
end

end

function [x, w, v] = lobpts_(n, dom)

if ( n == 1 )
    x = 0;
    w = 2;
    v = 1;
else
    [x, w, v] = lobpts(n);
end
x = (x+1)/2*diff(dom) + dom(1);
w = w/2*diff(dom);

end

function b = recursive(d, n, alpha, family)
% The barycentric d-simplex coordinates for a multiindex alpha with length n,
% based on a 1D node family.

xn = family(n+1);
b = zeros(d, 1);
if ( isempty(xn) )
    return
end

if ( d == 2 )
    b = xn([alpha(1) ; alpha(2)]);
    return
end

weight = 0;
for i = 1:d
    alpha_noti = alpha([1:i-1 i+1:end]);
    n_noti = n - (alpha(i)-1);
    w = xn(n_noti+1);
    br = recursive(d-1, n_noti, alpha_noti, family);
    b(1:i-1)   = b(1:i-1)   + w * br(1:i-1);
    b(i+1:end) = b(i+1:end) + w * br(i:end);
    weight = weight + w;
end

b = b / weight;

end

function y = from_unit(x, domain)

switch lower(domain)
    case 'unit'
        y = x;
    case 'biunit'
        y = 2*x-1;
    case 'barycentric'
        z = 1 - sum(x, 2);
        y = [x z];
    case 'equilateral'
        d = size(x, 2);
        y = x - 1/(d+1);               % Shift centroid to zero
        y = unit_to_equilateral(d, y);
        y = 2*y;                       % Scale edge length to 2
    otherwise
        error('Not implemented.')
end

end

function x = unit_to_equilateral(d, x)

if ( d > 1 )
    % Move the top vertex over the centroid
    x(:,1:d-1) = x(:,1:d-1) + x(:,d)/d;
    % Make the projection onto the lesser dimensions equilateral
    x(:,1:d-1) = unit_to_equilateral(d-1, x(:,1:d-1));
    % Scale the vertical dimension
    x(:,d) = x(:,d) * sqrt((d+1)/(2*d));
end

end

function y = to_unit(x, domain)

switch lower(domain)
    case 'unit'
        y = x;
    case 'biunit'
        y = (x+1)/2;
    case 'barycentric'
        y = x(:,1:end-1);
    case 'equilateral'
        d = size(x, 2);
        y = x / 2;                     % Scale edge length to 1
        y = equilateral_to_unit(d, y);
        y = y + 1/(d+1);               % Shift the first vertex to zero
    otherwise
        error('Not implemented.');
end

end

function x = equilateral_to_unit(d, x)

if ( d > 1 )
    % Scale the vertical dimension
    x(:,d) = x(:,d) / sqrt((d+1)/(2*d));
    % Make the projection onto the lesser dimensions right-angled
    x(:,1:d-1) = equilateral_to_unit(d-1, x(:,1:d-1));
    % Move the top vertex over first vertex
    x(:,1:d-1) = x(:,1:d-1) - x(:,d)/d;
end

end
