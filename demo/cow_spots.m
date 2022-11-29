delta_v = 1e-3;
delta_u = 0.516*delta_v;
alpha = 0.899;
beta = -0.91;
gamma = -0.899;
tau1 = 0.02; tau2 = 0.2;
dt = 0.1;

Nu = @(u,v) alpha*u.*(1-tau1*v.^2) + v.*(1-tau2*u);
Nv = @(u,v) beta*v.*(1+alpha/beta*tau1*u.*v) + u.*(gamma+tau2*v);

dom = surfacemesh.fromRhino('../models/cow.csv', 8);
dom = resample(dom, 12);

rng(1)
pdo = struct('lap', -dt*delta_u, 'b', 1);
Lu = surfaceop(dom, pdo);
build(Lu)
pdo = struct('lap', -dt*delta_v, 'b', 1);
Lv = surfaceop(dom, pdo);
build(Lv)

%% Initial conditions
bb = boundingbox(dom);
f = randnfun3(0.2, bb);
uinit = surfacefun(@(x,y,z) f(x,y,z), dom);
f = randnfun3(0.2, bb);
vinit = surfacefun(@(x,y,z) f(x,y,z), dom);

%% Simulation
close all
t = 0;
u = uinit;
v = vinit;

doplot = @(u) chain(@()plot(u), @()axis('off'), @()material('dull'), ...
    @()lighting('gouraud'), @()colormap('turbo'), @()colorbar, @()caxis('auto'), ...
    @() view(180, -85), @()camorbit(20, 0, 'data', 'y'), @()camorbit(10, 0, 'data', 'x'), ...
    @()camorbit(-3, 0, 'data', [0 0 1]), @()camlight('headlight'));

clf, doplot(u), drawnow, shg

kend = 2000;
for k = 1:kend
    tic
    Lu.rhs = u + dt*Nu(u,v);
    Lv.rhs = v + dt*Nv(u,v);
    u = solve(Lu);
    v = solve(Lv);
    toc
    t = t + dt;
    if ( mod(k, 10) == 0 )
        fprintf('k = %d\n', k);
        clf, doplot(u), drawnow, shg
    end
end
