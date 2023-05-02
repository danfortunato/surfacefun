function dom = twisted_torus(n, nu, nv)

if ( nargin < 2 )
    nu = 8;
end

if ( nargin < 3 )
    nv = nu;
end


xL = cell(nu*nv, 1); xR = cell(nu*nv, 1); xB = cell(nu*nv, 1); xF = cell(nu*nv, 1);
yL = cell(nu*nv, 1); yR = cell(nu*nv, 1); yB = cell(nu*nv, 1); yF = cell(nu*nv, 1);
zL = cell(nu*nv, 1); zR = cell(nu*nv, 1); zB = cell(nu*nv, 1); zF = cell(nu*nv, 1);

[nu, nv] = deal(nv, nu);
ubreaks = 0:nu;
vbreaks = linspace(0, 1, nv+1);
k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [xL{k}, yL{k}, zL{k}] = deal( 0+0*vv, vv, uu);
        [xR{k}, yR{k}, zR{k}] = deal( 1+0*vv, vv, uu);
        [xB{k}, yB{k}, zB{k}] = deal(vv,  0+0*vv, uu);
        [xF{k}, yF{k}, zF{k}] = deal(vv,  1+0*vv, uu);
        k = k+1;
    end
end

x = [xL ; xR ; xB ; xF];
y = [yL ; yR ; yB ; yF];
z = [zL ; zR ; zB ; zF];

for k = 1:length(x)
    [x{k}, y{k}, z{k}] = twist(x{k}, y{k}, z{k}, nu);
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = twist(x, y, z, nphi)

ns = 3;
theta0 = 0;
r = 0.3;
R = 1;

phi = 2*pi*z / nphi;
theta = theta0 + phi * ns/4;
u = sqrt(2) * (y - 0.5) * r;
v = sqrt(2) * (x - 0.5) * r;
x = (R + u.*cos(theta) + v.*sin(theta)) .* cos(phi);
y = (R + u.*cos(theta) + v.*sin(theta)) .* sin(phi);
z = v.*cos(theta) - u.*sin(theta);

end
