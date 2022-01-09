%% Geometry
n = 20;
%dom = surfacemesh.sphere(n, 2);
rng(1)
dom = surfacemesh.blob(n, 2);

%% Parameters
params = struct();

% Biological parameters:
params.delta = sqrt(0.1);
params.alpha = 1.5;
params.beta  = 0.1;
params.nu    = 100;
params.gamma = 0.23;

% Discretization parameters:
params.tend = 20;
params.dt   = 0.1;
params.scheme = 'bdf1';
params.useHeaviside     = false;
params.trueConservation = true;

% Visualization parameters:
params.quiet    = false;
params.movie    = true;
params.colormap = 'jet';
% params.movfile  = VideoWriter('prolate_polarization.mp4', 'MPEG-4');

%% Initial conditions:
init = 'random';

switch lower(init)
    case 'random'
        %rng(0)
        count = 0;
        i = 1;
        rng(count + 11522 +  200*i)
        %U = randnfun3(0.6, [-b b -b b -a a]);
        u = randnfun3(2, [-3 3 -3 3 -3 3]);
        u = u - min3(u) + 1e-8;  % Restrict to positive values
        u = 0.4*u ./ max3(u);        % Normalize
%         U = randnfunsphere(0.1); % Random spherical harmonic expansion
%         U = U - min2(U) + 1e-8;  % Restrict to positive values
%         U = U ./ max2(U) / 1.5;        % Normalize
    case 'gaussian'
        u = @(x,y,z) exp(-4*(x.^2+(y+0.7).^2+(z-0.6).^2));
    otherwise
        error('Unknown initial condition.');
end

u = surfacefun(@(x,y,z) u(x,y,z), dom);

%% Simulation
tic
v = SurfaceDiffusion(u, params);
toc
