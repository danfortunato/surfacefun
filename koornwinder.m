function [K, Kx, Ky] = koornwinder(n, x, y)
%KOORNWINDER   Evaluate Koornwinder polynomials.
%   K = KOORNWINDER(N, X, Y) constructs the Vandermonde matrix K mapping
%   NPOLY = (N+1)(N+2)/2 coefficients of an order-N Koornwinder polynomial
%   expansion to values of the expansion at the given points (X,Y). The points
%   (X,Y) must lie in the reference triangle with vertices (0,0), (1,0), (0,1).
%   If X and Y are M x 1 vectors, then K is an M x NPOLY matrix.
%
%   [K, KX, KY] = KOORNWINDER(N, X, Y) also constructs Vandermonde matrices KX
%   and KY for the x- and y-derivatives of the Koornwinder polynomials.
%
%   The Koornwinder polynomials K_{n,k}(x,y) are defined as
%
%      K_{n,k}(x,y) = P_{n-k}^(2k+1,0)(2x-1) * P_k(2y/(1-x)-1) * (1-x)^k
%
%   for 0<=n<=N, 0<=k<=n, and (x,y) in the reference triangle.

if ( nargin < 2 )
    [x, y] = trianglepts(n+1);
end

x = x(:);
y = y(:);

m = numel(x);
npoly = (n+1)*(n+2)/2;
K  = zeros(m, npoly);
Kx = zeros(m, npoly);
Ky = zeros(m, npoly);

% Compute the scaled Legendre polynomials P_k(2*(y/(1-x))-1) (1-x)^k
z = 2*y-(1-x);
leg  = zeros(m, n+1);
legx = zeros(m, n+1);
legy = zeros(m, n+1);
leg(:,1) = 1; legx(:,1) = 0; legy(:,1) = 0;
leg(:,2) = z; legx(:,2) = 1; legy(:,2) = 2;
for k = 2:n
    leg(:,k+1)  = ((2*k-1)*z.*leg(:,k) - (k-1)*(1-x).^2.*leg(:,k-1)) / k;
    legx(:,k+1) = ((2*k-1)*(  leg(:,k) + z.*legx(:,k)) - (k-1)*((1-x).^2.*legx(:,k-1) - 2*(1-x).*leg(:,k-1))) / k;
    legy(:,k+1) = ((2*k-1)*(2*leg(:,k) + z.*legy(:,k)) - (k-1)*(1-x).^2.*legy(:,k-1)) / k;
end

i = 1;
for nn = 0:n
    for kk = 0:nn
        scl = sqrt(1/(2*(2*kk+1)*(nn+1)));
        jac  =    jacobi(nn-kk, 2*kk+1, 0, 2*x-1);
        jacx = 2*djacobi(nn-kk, 2*kk+1, 0, 2*x-1);
        K(:,i)  = leg(:,kk+1) .*jac/scl;
        Ky(:,i) = legy(:,kk+1).*jac/scl;
        Kx(:,i) = legx(:,kk+1).*jac/scl + leg(:,kk+1).*jacx/scl;
        i = i+1;
    end
end

end

function djac = djacobi(n, a, b, x)

djac = zeros(size(x));
if ( n > 0 )
    djac = 0.5*(n+a+b+1)*jacobi(n-1,a+1,b+1,x);
end

end
