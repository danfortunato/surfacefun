function dom = sharptorus(n, nu, nv)

if ( nargin < 2 )
    nu = 8;
end

if ( nargin < 3 )
    nv = nu;
end

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalTorus(u, v)

theta = u;
phi = v;

x = phi/(2*pi);

y1 = 2*(2*x-1).*(2*x-1) - 1;
y2 = 15*x.*(x-0.6).*(x-1);

z = y1 * cos(3.75) - y2 * sin(3.75);
r = 2.5 + y1 * sin(3.75) + y2 * cos(3.75);

x = r .* cos(theta);
y = r .* sin(theta);

end
