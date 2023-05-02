function dom = mobius(n, nu, nv)

if ( nargin < 2 )
    nu = 8;
end

if ( nargin < 3 )
    nv = 2;
end

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(-1, 1, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalMobiusStrip(uu, vv);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);
%dom = surfacemesh([x ; x], [y ; y], [z ; z]);

end

function [x, y, z] = evalMobiusStrip(u, v)

r = 1;
k = 1;
x = (r+v/2.*cos(k*u/2)).*cos(u);
y = (r+v/2.*cos(k*u/2)).*sin(u);
z = v/2.*sin(k*u/2);

end
