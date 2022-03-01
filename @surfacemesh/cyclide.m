function dom = cyclide(n, nu, nv)

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
        [x{k}, y{k}, z{k}] = evalCyclide(uu, vv);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalCyclide(u, v)

a = 1;
b = 0.98;
c = 0.199;
d = 0.3;

a = 2;
b = 1.9;
c = sqrt(a^2-b^2);
d = 1;

x = (d*(c-a*cos(u).*cos(v)) + b^2*cos(u)) ./ (a-c*cos(u).*cos(v));
y = b*sin(u).*(a-d*cos(v)) ./ (a-c*cos(u).*cos(v));
z = b*sin(v).*(c*cos(u)-d) ./ (a-c*cos(u).*cos(v));

end
