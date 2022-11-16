function L = initialize(op, dom, rhs)
%INITIALIZE   Initialize an array of LEAF objects.
%   L = SURFACEOP.LEAF.INITIALIZE(OP, DOM) returns a cell array L of LEAF
%   objects which contain the solution and D2N operators for Poisson's
%   equation on the domain DOM with zero righthand side.
%
%   L = SURFACEOP.LEAF.INITIALIZE(OP, DOM, RHS) is as above, but with the
%   righthand side RHS, which may be a scalar or a function handle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( isempty(dom) )
    L = [];
    return
end

assert(isa(dom, 'surfacemesh'), 'Invalid domain.');

if ( nargin < 3 )
    % Default to homogeneous problem:
    rhs = 0;
end

numPatches = length(dom);
n = size(dom.x{1}, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% DEFINE REFERENCE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%

[X, Y] = chebpts2(n);             % Chebyshev points and grid.
ii = abs(X) < 1 & abs(Y) < 1;     % Interior indices.
ee = ~ii;                         % Boundary indices.
numBdyPts = sum(ee(:));
numIntPts = sum(ii(:));

% Skeleton mappings
nskel = n-2;
numSkelPts = 4*nskel;
S2L = skel2leaf(n, nskel);
L2S = leaf2skel(nskel, n);
S2L = sparse(S2L);
L2S = sparse(L2S);
xskel = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(n, 2);
B = barymat(xskel, xleaf, wleaf);
w = chebtech1.quadwts(nskel); w = w(:);
wskel = [w ; w ; w ; w];

% Skeleton indices for each side
leftSkel  = 1:nskel;
rightSkel = nskel+1:2*nskel;
downSkel  = 2*nskel+1:3*nskel;
upSkel    = 3*nskel+1:4*nskel;

% Compute binormal vectors on the skeleton
[NL, NR, ND, NU] = binormals(dom);
NL = pagemtimes(B, NL);
NR = pagemtimes(B, NR);
ND = pagemtimes(B, ND);
NU = pagemtimes(B, NU);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = diffmat(n);
I = eye(n);
II = kron(I, I);
Du = kron(D, I);
Dv = kron(I, D);
opfields = fieldnames(op).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT RHS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define scalar RHSs:
if ( isnumeric(rhs) && isscalar(rhs) )
    rhs = repmat(rhs, numIntPts, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% SOLVE LOCAL PROBLEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each patch needs a different solution operator.

% Initialize
L = cell(numPatches, 1);

% Loop over each patch:
for k = 1:numPatches

    % Define the left and right edges for this patch:
    x = dom.x{k};
    y = dom.y{k};
    z = dom.z{k};
    edges = [ x(1,1) y(1,1) z(1,1) x(n,1) y(n,1) z(n,1) nskel ;  % "Left" side
              x(1,n) y(1,n) z(1,n) x(n,n) y(n,n) z(n,n) nskel ;  % "Right" side
              x(1,1) y(1,1) z(1,1) x(1,n) y(1,n) z(1,n) nskel ;  % "Down" side
              x(n,1) y(n,1) z(n,1) x(n,n) y(n,n) z(n,n) nskel ]; % "Up" side

    % Evaluate non-constant RHSs if required:
    if ( isa(rhs, 'function_handle') )
        rhs_eval = feval(rhs, x(ii), y(ii), z(ii));
    elseif ( isa(rhs, 'surfacefun') )
        rhs_eval = rhs.vals{k}(ii);
    elseif ( iscell(rhs) )
        rhs_eval = rhs{k};
    else
        rhs_eval = rhs;
    end

    Dx = dom.ux{k}(:).*Du + dom.vx{k}(:).*Dv;
    Dy = dom.uy{k}(:).*Du + dom.vy{k}(:).*Dv;
    Dz = dom.uz{k}(:).*Du + dom.vz{k}(:).*Dv;
    Dxx = Dx^2;
    Dyy = Dy^2;
    Dzz = Dz^2;
    Dxy = Dx*Dy;
    Dyz = Dy*Dz;
    Dxz = Dx*Dz;
    if ( dom.singular(k) )
        J = dom.J{k}(:);
        Jx = Dx*J; Jy = Dy*J; Jz = Dz*J;
        Dxx = J.*Dxx - Jx.*Dx;
        Dyy = J.*Dyy - Jy.*Dy;
        Dzz = J.*Dzz - Jz.*Dz;
        Dxy = J.*Dxy - Jx.*Dy;
        Dyz = J.*Dyz - Jy.*Dz;
        Dxz = J.*Dxz - Jz.*Dx;
        Dx = J.^2.*Dx;
        Dy = J.^2.*Dy;
        Dz = J.^2.*Dz;
        II = J.^3.*II;
    end

    opk = op;
    for name = opfields
        name = name{1};
        if ( isa(opk.(name), 'function_handle') )
            opk.(name) = feval(opk.(name), x, y, z);
        elseif ( isa(opk.(name), 'surfacefun') )
            opk.(name) = opk.(name).vals{k};
        end
    end

    A = zeros(n^2);
    if ( opk.dxx ~= 0 ), A = A + opk.dxx(:).*Dxx; end
    if ( opk.dyy ~= 0 ), A = A + opk.dyy(:).*Dyy; end
    if ( opk.dzz ~= 0 ), A = A + opk.dzz(:).*Dzz; end
    if ( opk.dxy ~= 0 ), A = A + opk.dxy(:).*Dxy; end
    if ( opk.dyz ~= 0 ), A = A + opk.dyz(:).*Dyz; end
    if ( opk.dxz ~= 0 ), A = A + opk.dxz(:).*Dxz; end
    if ( opk.dx  ~= 0 ), A = A + opk.dx(:).*Dx;   end
    if ( opk.dy  ~= 0 ), A = A + opk.dy(:).*Dy;   end
    if ( opk.dz  ~= 0 ), A = A + opk.dz(:).*Dz;   end
    if ( opk.b   ~= 0 ), A = A + opk.b(:).*II;    end

    % Construct solution operator:
    if ( dom.singular(k) )
        dA = decomposition(A(ii,ii), 'cod');
        %dA = decomposition(A(ii,ii));
        Ainv = @(u) dA \ (dom.J{k}(ii).^3 .* u);
        %S = Ainv([-A(ii,ee), rhs_eval]);
        S = dA \ ([-A(ii,ee), dom.J{k}(ii).^3.*rhs_eval]);
    else
        dA = matlab.internal.decomposition.DenseLU(A(ii,ii));
        Ainv = @(u) solve(dA, u, false);
        S = Ainv([-A(ii,ee), rhs_eval]);
    end

    % Append boundary points to solution operator:
    tmpS = zeros(n^2, numBdyPts+1);
    tmpS(ii,:) = S;
    tmpS(ee,:) = eye(numBdyPts, numBdyPts+1);
    S = [tmpS(:,1:end-1) * S2L, tmpS(:,end)];

    % Construct normal derivative operator:
    nl = NL(:,:,k);
    nr = NR(:,:,k);
    nd = ND(:,:,k);
    nu = NU(:,:,k);
    dx = L2S * Dx(ee,:);
    dy = L2S * Dy(ee,:);
    dz = L2S * Dz(ee,:);
    normal_d = zeros(numSkelPts, n^2);
    normal_d(leftSkel,:)  = nl(:,1).*dx(leftSkel,:)  + nl(:,2).*dy(leftSkel,:)  + nl(:,3).*dz(leftSkel,:);
    normal_d(rightSkel,:) = nr(:,1).*dx(rightSkel,:) + nr(:,2).*dy(rightSkel,:) + nr(:,3).*dz(rightSkel,:);
    normal_d(downSkel,:)  = nd(:,1).*dx(downSkel,:)  + nd(:,2).*dy(downSkel,:)  + nd(:,3).*dz(downSkel,:);
    normal_d(upSkel,:)    = nu(:,1).*dx(upSkel,:)    + nu(:,2).*dy(upSkel,:)    + nu(:,3).*dz(upSkel,:);

    % Construct the D2N map:
    D2N = normal_d * S;

    % The D2N map needs to be scaled on each side (e.g. when being
    % merged) to account for the Jacobian scaling which has been
    % factored out of the coordinate derivative maps. This scaling
    % is not known until the merge stage, as it depends on the
    % scaling of the neighboring patch.
    if ( dom.singular(k) )
        D2N_scl = cell(4, 1);
        J = dom.J{k}.^3;
        Jss = L2S * J(ee);
        D2N_scl{1} = Jss(leftSkel);  % Left
        D2N_scl{2} = Jss(rightSkel); % Right
        D2N_scl{3} = Jss(downSkel);  % Down
        D2N_scl{4} = Jss(upSkel);    % Up
    else
        D2N_scl = {ones(nskel,1) ; ones(nskel,1) ; ones(nskel,1) ; ones(nskel,1)};
    end

    % Extract the particular solution to store separately:
    u_part = S(:,end); S = S(:,1:end-1);
    du_part = D2N(:,end); D2N = D2N(:,1:end-1);
    
    JJ = L2S * sqrt(dom.J{k}(ee));
    ww = wskel .* JJ;

    % Assemble the patch:
    xee = x(ee);
    yee = y(ee);
    zee = z(ee);
    xyz = L2S * [xee yee zee];
    L{k} = surfaceop.leaf(dom, k, S, D2N, D2N_scl, u_part, du_part, edges, xyz, ww, Ainv, normal_d);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nl, nr, nd, nu] = binormals(dom)
%BINORMALS   Compute the binormal vectors for a surfacemesh.

n = size(dom.x{1}, 1);

xu = cat(3, dom.xu{:}); xv = cat(3, dom.xv{:});
yu = cat(3, dom.yu{:}); yv = cat(3, dom.yv{:});
zu = cat(3, dom.zu{:}); zv = cat(3, dom.zv{:});

% Normal vectors to the surface (unnormalized)
nl = -[xu(:,1,:)   yu(:,1,:)   zu(:,1,:)];
nr =  [xu(:,n,:)   yu(:,n,:)   zu(:,n,:)];
nd = -[xv(1,:,:) ; yv(1,:,:) ; zv(1,:,:)]; nd = pagetranspose(nd);
nu =  [xv(n,:,:) ; yv(n,:,:) ; zv(n,:,:)]; nu = pagetranspose(nu);

% Tangent vectors to the element boundary (normalized)
tl = normalize([xv(:,1,:)   yv(:,1,:)   zv(:,1,:)]);
tr = normalize([xv(:,n,:)   yv(:,n,:)   zv(:,n,:)]);
td = normalize(pagetranspose([xu(1,:,:) ; yu(1,:,:) ; zu(1,:,:)]));
tu = normalize(pagetranspose([xu(n,:,:) ; yu(n,:,:) ; zu(n,:,:)]));

% Binormal vectors (normalized)
nl = normalize(nl - tl .* sum(nl.*tl, 2));
nr = normalize(nr - tr .* sum(nr.*tr, 2));
nd = normalize(nd - td .* sum(nd.*td, 2));
nu = normalize(nu - tu .* sum(nu.*tu, 2));

end

function v = normalize(v)

v = v ./ sqrt(v(:,1,:).^2 + v(:,2,:).^2 + v(:,3,:).^2);

end

function P = skel2leaf(nleaf, nskel)
%SKEL2LEAF   Boundary interpolation matrix.
%   SKEL2LEAF(NLEAF, NSKEL) returns the (4*NLEAF-4) x 4*NSKEL matrix that
%   maps 4 pieces of length-NSKEL first-kind boundary values to 4*NLEAF-4
%   second-kind boundary values, including the corners. At each corner, the
%   average of the two interpolated values is used.

[xskel, ~, wskel] = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(nleaf, 2);
B = barymat(xleaf, xskel, wskel);

% Skeleton indices for each side
leftSkel  = 1:nskel;
rightSkel = nskel+1:2*nskel;
downSkel  = 2*nskel+1:3*nskel;
upSkel    = 3*nskel+1:4*nskel;

% Leaf indices for each side
leftLeaf  = 1:nleaf;
rightLeaf = 3*nleaf-3:4*nleaf-4;
upLeaf    = [nleaf:2:3*nleaf-4 4*nleaf-4];
downLeaf  = [1 nleaf+1:2:3*nleaf-3];

P = zeros(4*nleaf-4, 4*nskel);
P(leftLeaf,  leftSkel)  = B;
P(rightLeaf, rightSkel) = B;
P(downLeaf,  downSkel)  = B;
P(upLeaf,    upSkel)    = B;

% Average the corners:
corners = [1 nleaf 3*nleaf-3 4*nleaf-4];
P(corners,:) = P(corners,:)/2;

end

function P = leaf2skel(nskel, nleaf)
%LEAF2SKEL   Boundary interpolation matrix.
%   LEAF2SKEL(NSKEL, NLEAF) returns the 4*NSKEL x (4*NLEAF-4) matrix that
%   maps 4*NLEAF-4 second-kind boundary values to 4 pieces of length-NSKEL
%   first-kind boundary values.

[xskel, ~, wskel] = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(nleaf, 2);
B = barymat(xskel, xleaf, wleaf);

% Skeleton indices for each side
leftSkel  = 1:nskel;
rightSkel = nskel+1:2*nskel;
downSkel  = 2*nskel+1:3*nskel;
upSkel    = 3*nskel+1:4*nskel;

% Leaf indices for each side
leftLeaf  = 1:nleaf;
rightLeaf = 3*nleaf-3:4*nleaf-4;
upLeaf    = [nleaf:2:3*nleaf-4 4*nleaf-4];
downLeaf  = [1 nleaf+1:2:3*nleaf-3];

P = zeros(4*nskel, 4*nleaf-4);
P(leftSkel,  leftLeaf)  = B;
P(rightSkel, rightLeaf) = B;
P(downSkel,  downLeaf)  = B;
P(upSkel,    upLeaf)    = B;

end
