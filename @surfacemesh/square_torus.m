function dom = square_torus(n, nu, nv)

if ( nargin < 2 )
    nu = 8;
end

if ( nargin < 3 )
    nv = nu;
end

router = 1;
rinner = 0.5;
zb = -0.25;
zt = 0.25;

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);

%% Outer ring
ubreaks = linspace(zb, zt, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);
k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv, router);
        k = k+1;
    end
end
xouter = x;
youter = y;
zouter = z;


%% Inner ring
k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv, rinner);
        k = k+1;
    end
end
xinner = x;
yinner = y;
zinner = z;

%% Top
ubreaks = linspace(rinner, router, nu+1);
k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalAnnulus(uu, vv, zt);
        k = k+1;
    end
end
xtop = x;
ytop = y;
ztop = z;

%% Bottom
k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalAnnulus(uu, vv, zb);
        k = k+1;
    end
end
xbottom = x;
ybottom = y;
zbottom = z;

x = [xouter ; xtop ; xinner ; xbottom];
y = [youter ; ytop ; yinner ; ybottom];
z = [zouter ; ztop ; zinner ; zbottom];

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalTorus(u, v, r)

phi = v;
x = r .* cos(phi);
y = r .* sin(phi);
z = u;

end

function [x, y, z] = evalAnnulus(u, v, z)

r = u;
theta = v;
x = r .* cos(theta);
y = r .* sin(theta);
z = z + 0*u;

end
