function [u, dt] = SurfaceDiffusion(u, params)
%SORDIFFUSION   Solve the diffusion problem on a surface of revolution.
%   U = SORDIFFUSION(U, PARAMS) solves the PDE
%
%      u_t = delta^2*lap_s(u) - u
%            + (beta + u^nu/(u^nu + gamma^nu))*(1 - 2*pi*alpha*mean2(u))
%
%   on a surface of revolution. Here, lap_s is the Laplace-Beltrami
%   operator on the surface and mean2(u) = (int_{surface} u dS) /
%   (int_{surface} dS). U should be given as a matrix of double Fourier
%   coefficients.
%
%   The parameters are:
%
%      Biological parameters
%      ----------------------
%      PARAMS.DELTA:            Diffusivity
%      PARAMS.ALPHA:            Coupling parameter
%      PARAMS.BETA:             Constitutive rate
%      PARAMS.GAMMA:            Membrane recruitment threshold
%      PARAMS.NU:               Cooperativity parameter
%
%      Discretization parameters
%      -------------------------
%      PARAMS.DT:               Time step
%      PARAMS.TEND:             End time
%      PARAMS.SCHEME:           Time-stepping scheme ('bdf1', 'bdf2', or 'bdf4')
%      PARAMS.USEHEAVISIDE:     Replace the nonlinearity with a Heaviside
%      PARAMS.TRUECONSERVATION: Enforce conservation according to the true volume
%
%      Visualization parameters
%      ------------------------
%      PARAMS.QUIET:            Flag to keep quiet
%      PARAMS.MOVIE:            Flag to plot a movie during simulation
%      PARAMS.COLORMAP:         Colormap for plotting
%      PARAMS.MOVFILE:          Filename of movie to output

% Biological parameters:
delta = params.delta;
alpha = params.alpha;
beta  = params.beta;
nu    = params.nu;
gamma = params.gamma;
gamma_nu = gamma^nu;

% Discretization parameters:
tend  = params.tend;
dt    = params.dt;

% Make sure dt divides into tend nicely
nsteps = ceil(tend/dt);
dt = tend/nsteps;

switch ( lower(params.scheme) )
    case 'bdf1'
        cscl = -1/(dt*delta^2);
        step = @bdf1;
        multistep = 1;
    case 'bdf2'
        cscl = -3/(2*dt*delta^2);
        step = @bdf2;
        multistep = 2;
    case 'bdf4'
        cscl = -25/(12*dt*delta^2);
        step = @bdf4;
        multistep = 4;
    otherwise
        error('Unknown time-stepping scheme.');
end

% Build the Laplace-Beltrami solver
pdo = struct('lap', 1, 'b', cscl);
L = surfaceop(u.domain, pdo);
build(L)

outputMovie = isfield(params, 'movfile') && ~isempty(params.movfile);
if ( outputMovie )
    open(params.movfile);
end

if ( params.movie )
    % Plot the initial condition:
    clf
    plot(u), hold on, plot(u.domain)
    colormap(params.colormap)
    setupfigure(0)
    shg

    if ( outputMovie )
        frame = getframe(gcf);
        writeVideo(params.movfile, frame);
    end

    if ( ~params.quiet )
        fprintf('   Press any key when ready.\n\n')
        pause
    end
end

% The true conservation law changes depending on the volume of the
% geometry, which has the effect of rescaling alpha.
vol = volume(u.domain);
if ( params.trueConservation )
    % Conservation based on the true volume
    fac = 3/2 * vol;
else
    % Conservation based on the volume of a sphere
    fac = 2*pi;
end

    function u = N_direct(u)
        u_nu = u.^nu;
        u = (1-alpha/fac*integral2(u)) * (beta+u_nu./(u_nu+gamma_nu)) - u;
    end

    function u = N_heaviside(u)
        v = util.coeffs2valsDbl(u);
        mu = params.nlat/(5*log(params.nlat));
        v = beta + 1./(1+exp(-2*mu*(v-gamma)));
        U2 = util.vals2coeffsDbl(v);
        scl = 1 - alpha/fac * L.integral2(u);
        u = scl*U2 - u;
    end

if ( params.useHeaviside )
    N = @N_heaviside;
else
    N = @N_direct;
end

t = 0;
u_old = cell(multistep-1, 1);
if ( multistep > 1 )
    % Start the multistep method with a few steps of IMEX BDF1
    cscl1 = -1/(dt*delta^2);
    pdo = struct('lap', 1, 'b', cscl1);
    L1 = surfaceop(u.domain, pdo);
    build(L1)
    for k = 1:multistep-1
        rhs = bdf1(N, dt, u);
        u_old{multistep-k} = u;
        L1.rhs = cscl1 * rhs;
        u = solve(L1);
        t = t + dt;
        nsteps = nsteps-1;
    end
end

for k = 1:nsteps

    [rhs, u_old] = step(N, dt, u, u_old);
    L.rhs = cscl * rhs;
    u = solve(L);
    t = t + dt;

    if ( params.movie )
        clf
        plot(u), hold on, plot(u.domain)
        setupfigure(t)
        drawnow

        if ( outputMovie )
            frame = getframe(gcf);
            writeVideo(params.movfile, frame);
        end
    end
end

if ( outputMovie )
    close(params.movfile);
end

end

%% Time-stepping schemes

function [rhs, u_old] = bdf1(N, dt, u, u_old)
    rhs = u + dt*N(u);
end

function [rhs, u_old] = bdf2(N, dt, u, u_old)
    rhs = 4/3*u - 1/3*u_old{1} + dt*(4/3*N(u) - 2/3*N(u_old{1}));
    u_old{1} = u;
end

function [rhs, u_old] = bdf4(N, dt, u, u_old)
    rhs = 48/25*u - 36/25*u_old{1} + 16/25*u_old{2} - 3/25*u_old{3} + ...
          dt*(48/25*N(u) - 72/25*N(u_old{1}) + 48/25*N(u_old{2}) - 12/25*N(u_old{3}));
    u_old(2:3) = u_old(1:2);
    u_old{1} = u;
end

%% Plotting utilities

function setupfigure(t)
    axis off
    title(sprintf('$t = %.2f$', t), 'Interpreter', 'Latex', 'FontSize', 18)
    colorbar('FontSize', 14)
end
