function dom = torus(n, nu, nv)

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

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.45 0];

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));

for i = -1:2
    for j = -1:2
        ii = i+2;
        jj = j+2;
        x = x + d(ii,jj)*cos(v).*cos((1-i)*u+j*v);
        y = y + d(ii,jj)*sin(v).*cos((1-i)*u+j*v);
        z = z + d(ii,jj)*sin((1-i)*u+j*v);
    end
end
 
end
