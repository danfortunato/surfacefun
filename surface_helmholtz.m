rng(0)
n = 16;             % Number of nodes per dimension for each mesh element
nref = 3;           % Number of levels of uniform mesh refinement
rect = [-1 1 -1 1]; % Bounding box of compactly supported region
zk = 30;            % Helmholtz wavenumber
eta = zk;           % Impedance parameter

% Incoming plane wave
uincf    = @(x,y,z) exp(1i*zk*x);
uincf_dx = @(x,y,z) 1i*zk*exp(1i*zk*x);
uincf_dy = @(x,y,z) 0*x;

%% Make a compactly supported perturbation of a flat surface
% Build a uniformly refined flat square
dom = surfacemesh.square(n, nref, rect);
x = dom.x;
y = dom.y;
z = dom.z;

% Perturb it in the z-direction according to a compactly supported function
% f = randnfun2(0.2, rect); f = f ./ max(abs(min2(f)), max2(f));
% f = @(x,y) cos(4*pi*x).*cos(4*pi*y);
f = @(x,y) 0.5+0*x;
b = @(x,y,z) f(x,y).*exp(-50*(x.^2+y.^2));
for k = 1:length(dom)
    z{k} = b(x{k},y{k},z{k});
end
dom = surfacemesh(x, y, z);

%% Construct the HPS solver for the total field and get the interior DtN operator
pdo = [];
pdo.lap = 1;
pdo.c = zk^2;
L = surfaceop(dom, pdo, method='ItI', eta=eta);
DtN_int = L.DtN();
xyz = L.patches{1}.xyz;

%% Convert from HPS panels to chunkie panels
% Collect indices corresponding to each side of the outer boundary
leftIdx  = find(abs(xyz(:,1) - rect(1)) < 1e-10);
rightIdx = find(abs(xyz(:,1) - rect(2)) < 1e-10);
downIdx  = find(abs(xyz(:,2) - rect(3)) < 1e-10);
upIdx    = find(abs(xyz(:,2) - rect(4)) < 1e-10);

% Reorder the HPS skeleton to match chunkie's panel ordering
nskel = n-2;
nside = nskel * 2^nref;
P = speye(4*nside);
P(:, [downIdx; rightIdx; flip(upIdx); flip(leftIdx)]) = P;

% Interpolate from Chebyshev panels (HPS) to Gauss panels (chunkie)
[xc, ~, vc] = chebpts(nskel, 1);
[xl, ~, vl] = legpts(nskel);
C2L = barymat(xl, xc, vc); C2L = repmat({C2L}, 4*2^nref, 1); C2L = matlab.internal.math.blkdiag(C2L{:});
L2C = barymat(xc, xl, vl); L2C = repmat({L2C}, 4*2^nref, 1); L2C = matlab.internal.math.blkdiag(L2C{:});
Q = L2C * P';
P = P * C2L;

DtN_int = P * DtN_int * Q;
xyz     = P * xyz;

%% Create the boundary as a chunker and construct layer potentials
chnkr = squarechunker(nskel, nref, rect);
Skern = kernel('helmholtz', 's', zk);
Dkern = kernel('helmholtz', 'd', zk);
S = chunkermat(chnkr, Skern);
D = chunkermat(chnkr, Dkern);
I = eye(chnkr.npt);

%% Solve for the scattered field on the boundary
% Evaluate the incoming field and its normal derivative
uinc    = uincf(xyz(:,1), xyz(:,2), xyz(:,3));
uinc_dx = uincf_dx(xyz(:,1), xyz(:,2), xyz(:,3));
uinc_dy = uincf_dy(xyz(:,1), xyz(:,2), xyz(:,3));
uinc_dn = [ -uinc_dy(1:nside)           ;  % Bottom
             uinc_dx(nside+1:2*nside)   ;  % Right
             uinc_dy(2*nside+1:3*nside) ;  % Top
            -uinc_dx(3*nside+1:4*nside) ]; % Left

% Solve for the scattered field and compute its normal derivative
A = I/2 - D + S*DtN_int;
rhs = S * (uinc_dn - DtN_int*uinc);
uscat = A \ rhs;
uscat_dn = DtN_int * (uinc + uscat) - uinc_dn;

% Construct impedance boundary data for the total field
dir = uinc    + uscat;
neu = uinc_dn + uscat_dn;
imp = Q * (neu + 1i*eta*dir);
u = L.solve(imp);

%% Evaluate outside the compact region
m = 100;
padx = 0.2*diff(rect(1:2));
pady = 0.2*diff(rect(3:4));
[xx, yy] = meshgrid(linspace(rect(1)-padx, rect(2)+padx, m), ...
                    linspace(rect(3)-pady, rect(4)+pady, m));
ii = xx < rect(1) | xx > rect(2) | yy < rect(3) | yy > rect(4);
targ = [xx(ii).' ; yy(ii).'];
opts = [];
opts.accel = true;
D_uscat  = chunkerkerneval(chnkr, Dkern, uscat,    targ, opts);
S_duscat = chunkerkerneval(chnkr, Skern, uscat_dn, targ, opts);
vals = D_uscat - S_duscat;

vv = nan(m);
vv(ii) = vals + uincf(xx(ii),yy(ii));

%% Plot
figure
plot(real(u))
hold on
pcolor(xx, yy, real(vv))
view(3)
shading interp
colormap turbo
camlight
colorbar
alignfigs
