% Ginzburg-Landau equation on a surface:
%    u_t = d (1+bi) lap(u) + u - (1+ci) u |u|^2

b = -0.3;
c = 1.25;
d = 1e-3;
delta = d*(1+b*1i);
dt = 0.2;
N = @(u) u - (1+c*1i)*u.*(abs(u).^2);

n = 16;
%dom = surfacemesh.sphere(n, 2);
dom = surfacemesh.blob(n, 2);
pdo = struct('lap', -dt*delta, 'b', 1);
L = surfaceop(dom, pdo);
build(L)

%%
%uinit = randnfunsphere(0.1);
uinit = randnfun3(0.5, );
uinit = uinit / norm(uinit, inf);
uinit = surfacefun(@(x,y,z) uinit(x,y,z), dom);

%%
close all
u = uinit;

plot(u)
colorbar
shg
pause

t = 0;
for k = 1:1000
    L.rhs = u + dt*N(u);
    u = solve(L);
    t = t + dt;
    if ( mod(k, 5) == 0 )
        k
        plot(u)
        colorbar
        drawnow
        shg
    end
end