function [u, v, w, curlfree, divfree] = hodge(f)
%HODGE   Hodge decomposition of a SURFACEFUNV.
%   [U, V, W] = HODGE(F) computes the Hodge decomposition of the
%   SURFACEFUNV F over the surface. This decomposes the vector field F into
%   scalar functions U and V and a harmonic vector field W such that
%
%      F = grad(U) + cross(N, grad(V)) + W,
%
%   where N = normal(domain(F)). The harmonic vector field W satisfies
%   div(W) = 0 and div(cross(N, W)) = 0. U and V are returned as
%   SURFACEFUNs and W is returned as a SURFACEFUNV.
%
%   [U, V, W, CURLFREE, DIVFREE] = HODGE(F) also returns the vector fields
%   CURLFREE = grad(U) and DIVFREE = cross(N, grad(V)) as SURFACEFUNVs.

pdo = [];
pdo.lap = 1;
dom = domain(f);
n = normal(dom);

% Curl-free component: Solve lap(u) = div(f)
L = surfaceop(dom, pdo, div(f));
u = solve(L);

% Divergence-free component: Solve lap(v) = -div(n x f)
L.rhs = -div(cross(n, f));
v = solve(L);

% Harmonic component: w = f - grad(u) - n x grad(v)
curlfree = grad(u);
divfree = cross(n, grad(v));
w = f - curlfree - divfree;

end
