% Compute the Hodge decomposition of a vector field f:
%
%   f = grad(u) + n x grad(v) + w

n = 20;
dom = surfacemesh.stellarator(n, 8, 24);

x0 = 0.2; y0 = 0.2; z0 = 0.2;
denom = @(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(3/2);
V = cross([0 1 1], surfacefunv(@(x,y,z) (x-x0)./denom(x,y,z), ...
                               @(x,y,z) (y-y0)./denom(x,y,z), ...
                               @(x,y,z) (z-z0)./denom(x,y,z), dom));

vn = normal(dom);
f = -cross(vn, vn, V);

%%
tic
[u, v, w, curlfree, divfree] = hodge(f);
toc

%%
close all
doplot = @(f,k,t) chain(@() figure(k), ...
                        @() plot(norm(f)), ...
                        @() hold('on'), ...
                        @() quiver(normalize(f), 0.6, 6, 'Color', 'k'), ...
                        @() axis('equal', 'off'), ...
                        @() colormap('turbo'), ...
                        @() colorbar, ...
                        @() title(t, 'FontSize', 16));
doplot(f, 1, 'Vector field'), hold on, plot(dom)
doplot(curlfree, 2, 'Curl-free component')
doplot(divfree, 3, 'Divergence-free component')
doplot(w, 4, 'Harmonic component')
alignfigs

norm(div(divfree))
