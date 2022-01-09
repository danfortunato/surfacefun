% Compute the Hodge decomposition of a vector field f:
%
%   f = grad(u) + n x grad(v) + w

n = 16;
dom = surfacemesh.torus(n, 8, 24);

x0 = 0.2; y0 = 0.2; z0 = 0.2;
denom = @(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(3/2);
V = cross([0 1 1], surfacefunv(@(x,y,z) (x-x0)./denom(x,y,z), ...
                               @(x,y,z) (y-y0)./denom(x,y,z), ...
                               @(x,y,z) (z-z0)./denom(x,y,z), dom));
vn = normal(dom);
f = -cross(vn, vn, V);

[u, v, w, curlfree, divfree] = hodge(f);

%%
close all
doplot = @(f,k) chain(@() figure(k), ...
                      @() plot(norm(f)), ...
                      @() hold('on'), ...
                      @() quiver(normalize(f), 0.3, 5), ...
                      @() axis('equal'), ...
                      @() colorbar);
doplot(f, 1)
doplot(curlfree, 2)
doplot(divfree, 3)
doplot(w, 4)
alignfigs
