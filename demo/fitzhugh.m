% Fitzhugh-Nagumo on a surface
%
% u_t = delta_u lap(u) + 1/alpha u (1 - u) (u - (v+b)/a);
% v_t = delta_v lap(v) + u - v

n = 16;
dom = surfacemesh.cyclide(n, 32, 16);

a = 0.75;
b = 0.02;
alpha = 0.02;
delta_u = 0.03;
delta_v = 0;
dt = 0.02;
Nu = @(u,v) 1/alpha*u.*(1-u).*(u-(v+b)/a);
Nv = @(u,v) u-v;

pdo = struct('lap', -dt*delta_u, 'b', 1);
Lu = surfaceop(dom, pdo);
build(Lu)
pdo = struct('lap', -dt*delta_v, 'b', 1);
Lv = surfaceop(dom, pdo);
build(Lv)

%% Initial conditions
uinit = surfacefun(@(x,y,z) 0.5*(1+tanh(2*x+y)), dom);
vinit = surfacefun(@(x,y,z) 0.5*(1-tanh(3*z)), dom);

%% Simulation
close all
u = uinit;
v = vinit;

doplot = @(u) chain(@()plot(u), @()view(60,40), @()material('dull'), ...
    @()lighting('gouraud'), @()camlight('headlight'), @()colorbar, @()caxis('auto'));

doplot(u), shg
t = 0;
for k = 1:10000
    Lu.rhs = u + dt*Nu(u,v);
    Lv.rhs = v + dt*Nv(u,v);
    u = solve(Lu);
    v = solve(Lv);
    t = t + dt;
    if ( mod(k,10) == 0 )
        k
        clf, doplot(u), drawnow, shg
    end
end