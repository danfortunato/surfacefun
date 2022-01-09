% Laplace-Beltrami

n = 20;
dom = surfacemesh.sphere(n, 2);

l = 2; m = 2;
sol = spherefun.sphharm(l, m);
sol = surfacefun(@(x,y,z) sol(x,y,z), dom);
f = -l*(l+1)*sol;

pdo = [];
pdo.lap = 1;
tic
L = surfaceop(dom, pdo, f);
toc

tic
build(L)
toc

tic
L.rhs = f;
toc

tic
u = solve(L);
toc

clf
plot(u - sol)
hold on
wireframe(dom)
axis equal
colorbar
shg
