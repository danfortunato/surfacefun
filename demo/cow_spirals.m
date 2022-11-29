delta = 5e-4;
c = 1.5;
dt = 0.03;
N = @(u) u - (1+c*1i)*u.*(abs(u).^2);

dom = surfacemesh.fromRhino('../models/cow.csv', 8);
dom = resample(dom, 12);

rng(1)
pdo = struct('lap', -dt*delta, 'b', 1);
L = surfaceop(dom, pdo);
L.build()

%% Initial conditions
bb = boundingbox(dom);
f = randnfun3(0.2, bb);
uinit = surfacefun(@(x,y,z) f(x,y,z), dom);

%% Simulation
close all
t = 0;
u = uinit;

doplot = @(u) chain(@()plot(u), @()axis('off'), @()material('dull'), ...
    @()lighting('gouraud'), @()colormap('turbo'), @()colorbar, @()caxis('auto'), ...
    @() view(180, -85), @()camorbit(20, 0, 'data', 'y'), @()camorbit(10, 0, 'data', 'x'), ...
    @()camorbit(-3, 0, 'data', [0 0 1]), @()camlight('headlight'));

clf, doplot(u), drawnow, shg

kend = 2000;
for k = 1:kend
    tic
    L.rhs = u + dt*N(u);
    u = solve(L);
    toc
    t = t + dt;
    if ( mod(k, 10) == 0 )
        fprintf('k = %d\n', k);
        clf, doplot(u), drawnow, shg
    end
end
