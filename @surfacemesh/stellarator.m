function dom = stellarator(n, nu, nv)

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
        [x{k}, y{k}, z{k}] = evalStellarator(uu, vv);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalStellarator(u, v)

Q = 3;
R = 1;
r0 = 0;
z0 = 0;
d = [0.15  0.09  0.00  0.00  0.00;
     0.00  1.00  0.03 -0.01  0.00;
     0.08  4.00 -0.01 -0.02  0.00;
     0.01 -0.28 -0.28  0.03  0.02;
     0.00  0.09 -0.03  0.06  0.00;
     0.01 -0.02  0.02  0.00 -0.02];

rz = zeros(size(u));
for j = -1:4
    for k = -1:3
        jj = j+2;
        kk = k+2;
        rz = rz + d(jj,kk)*exp(-1i*j*u+1i*k*Q*v);
    end
end
rz = exp(1i*u).*rz;
r1 = real(rz);
z1 = imag(rz);
r = r0 + R*(r1-r0);
z = z0 + R*(z1-z0);

[x, y] = pol2cart(v, r);

end
