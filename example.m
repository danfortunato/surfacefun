n = 20;
[x, y, z] = cubedsphere(n, 2);

idx = 1:length(x);
x = x(idx);
y = y(idx);
z = z(idx);

dom = struct();
K = length(x);
for k = 1:K
    dom(k).x = x{k};
    dom(k).y = y{k};
    dom(k).z = z{k};
end

l = 7; m = 3;
sol = spherefun.sphharm(l, m);
rhs = @(x,y,z) -l*(l+1)*sol(x,y,z);
S = SurfaceSEM(dom, rhs);
build(S)
xyz = S.patches{1}.xyz;
bc = sol(xyz(:,1), xyz(:,2), xyz(:,3));
u = S.solve(bc);

for k = 1:K
    %surf(x{k}, y{k}, z{k}, u{k})
    %surf(x{k}, y{k}, z{k}, sol(x{k}, y{k}, z{k}))
    surf(x{k}, y{k}, z{k}, u{k} - sol(x{k}, y{k}, z{k}))
    hold on
end
%shading interp
axis equal
colorbar
shg

%% Test with spherefun
rng(0)
f = randnfunsphere;
fx = diff(f, 1);
fy = diff(f, 2);
fz = diff(f, 3);
lapf = lap(f);

n = 20;
[x, y, z] = cubedsphere(n, 1);
idx = 11;
xx = x{idx};
yy = y{idx};
zz = z{idx};
ff = f(xx,yy,zz);

D = diffmat(n);
xu = xx * D.'; xv = D * xx;
yu = yy * D.'; yv = D * yy;
zu = zz * D.'; zv = D * zz;
E = xu.*xu + yu.*yu + zu.*zu;
G = xv.*xv + yv.*yv + zv.*zv;
F = xu.*xv + yu.*yv + zu.*zv;
J = E.*G - F.^2;
ux = (G.*xu-F.*xv)./J; vx = (E.*xv-F.*xu)./J;
uy = (G.*yu-F.*yv)./J; vy = (E.*yv-F.*yu)./J;
uz = (G.*zu-F.*zv)./J; vz = (E.*zv-F.*zu)./J;
I = speye(n);
Du = kron(D, I);
Dv = kron(I, D);
Dx = ux(:).*Du + vx(:).*Dv;
Dy = uy(:).*Du + vy(:).*Dv;
Dz = uz(:).*Du + vz(:).*Dv;

dx = reshape(Dx*ff(:), n, n);
dy = reshape(Dy*ff(:), n, n);
dz = reshape(Dz*ff(:), n, n);
dxx = reshape(Dx*Dx*ff(:), n, n);
dyy = reshape(Dy*Dy*ff(:), n, n);
dzz = reshape(Dz*Dz*ff(:), n, n);

%surf(xx, yy, zz, dx - fx(xx,yy,zz))
%surf(xx, yy, zz, dy - fy(xx,yy,zz))
%surf(xx, yy, zz, dz - fz(xx,yy,zz))
surf(xx, yy, zz, dxx + dyy + dzz - lapf(xx,yy,zz))
shading interp
axis equal
colorbar
shg


normalize = @(v) v ./ sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
nl = -normalize([xu(:,1)   yu(:,1)   zu(:,1)]);
nr =  normalize([xu(:,n)   yu(:,n)   zu(:,n)]);
nd = -normalize([xv(1,:).' yv(1,:).' zv(1,:).']);
nu =  normalize([xv(n,:).' yv(n,:).' zv(n,:).']);

tangent = normalize([xv(:,1) yv(:,1) zv(:,1)]);
nl = nl - tangent .* dot(nl, tangent, 2);
tangent = normalize([xv(:,n) yv(:,n) zv(:,n)]);
nr = nr - tangent .* dot(nr, tangent, 2);
tangent = normalize([xu(1,:).' yu(1,:).' zu(1,:).']);
nd = nd - tangent .* dot(nd, tangent, 2);
tangent = normalize([xu(n,:).' yu(n,:).' zu(n,:).']);
nu = nu - tangent .* dot(nu, tangent, 2);

nl = normalize(nl);
nr = normalize(nr);
nd = normalize(nd);
nu = normalize(nu);

hold on
quiver3(xx(:,1), yy(:,1), zz(:,1), nl(:,1), nl(:,2), nl(:,3), 'b')
quiver3(xx(:,n), yy(:,n), zz(:,n), nr(:,1), nr(:,2), nr(:,3), 'r')
quiver3(xx(1,:).', yy(1,:).', zz(1,:).', nd(:,1), nd(:,2), nd(:,3), 'g')
quiver3(xx(n,:).', yy(n,:).', zz(n,:).', nu(:,1), nu(:,2), nu(:,3), 'm')