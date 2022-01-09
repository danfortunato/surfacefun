% Gray-Scott equation on a surface:
%    u_t = e1 lap(u) + u - u^3

delta = 0.01;
dt = 0.1;
N = @(u) u - u.^3;

n = 16;
dom = surfacemesh.sphere(n, 2);
pdo = struct('lap', -dt*delta, 'b', 1);
L = surfaceop(dom, pdo);
build(L)

%%
close all
u = surfacefun(@(x,y,z) cos(cosh(5*x.*z)-10*y), dom);

plot(u)
colorbar
shg
pause

t = 0;
for k = 1:100
    L.rhs = u + dt*N(u);
    u = solve(L);
    t = t + dt;
    if ( mod(k,10) == 0 )
        k
        plot(u)
        colorbar
        drawnow
        shg
    end
end