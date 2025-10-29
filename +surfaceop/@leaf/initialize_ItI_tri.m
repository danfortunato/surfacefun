function L = initialize_ItI_tri(op, dom, eta, rhs)
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
npts = length(dom.x{1});
n = order(dom)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% DEFINE REFERENCE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%

[xx, yy] = trianglepts(n);
ee_left  = 1:n;
ee_down  = cumsum([1 (n:-1:2)]);
ee_hypot = cumsum([n (n-1:-1:1)]);
eeIdx = unique([ee_left ee_down ee_hypot]);
ee = false(size(xx));
ee(eeIdx) = true;    % Boundary indices
ii = ~ee;            % Interior indices
numBdyPts = sum(ee(:));
numIntPts = sum(ii(:));

% Impose the PDE on the interior nodes
ipde = ii;

% Impose the BCs on the interior of the boundary nodes
% TODO
ibc  = false(n, 1);
ibc(2:n-1) = true;

% Skeleton mappings
nskel = n-2;
numSkelPts = 3*nskel;
S2L = skel2leaf(n, nskel); % Don't sparsify for speed
L2S = leaf2skel(nskel, n);
B2S = bdy2skel(nskel, n); % TODO
xskel = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(n, 2);
B = barymat(xskel, xleaf, wleaf);
w = chebtech1.quadwts(nskel); w = w(:);
wskel = [0.5*w ; 0.5*w ; sqrt(2)/2*w];

% Skeleton indices for each side
leftSkel  = 1:nskel;
downSkel  = nskel+1:2*nskel;
hypotSkel = 2*nskel+1:3*nskel;

% Compute binormal vectors
[NL, ND, NH] = binormals(dom);

ux = reshape([dom.ux{:}], [npts numPatches]); vx = reshape([dom.vx{:}], [npts numPatches]);
uy = reshape([dom.uy{:}], [npts numPatches]); vy = reshape([dom.vy{:}], [npts numPatches]);
uz = reshape([dom.uz{:}], [npts numPatches]); vz = reshape([dom.vz{:}], [npts numPatches]);

ItI_scl0 = {ones(nskel,1) ; ones(nskel,1) ; ones(nskel,1)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[K, Ku, Kv] = koornwinder(n-1, xx, yy);
Du = Ku / K;
Dv = Kv / K;
II = eye(npts);

X = reshape([dom.x{:}], [npts numPatches]);
Y = reshape([dom.y{:}], [npts numPatches]);
Z = reshape([dom.z{:}], [npts numPatches]);

flags = structfun(@(f) ~(isscalar(f) && isnumeric(f) && f==0), op, 'UniformOutput', false);

for name = fieldnames(op).'
    name = name{1};
    if ( isa(op.(name), 'function_handle') )
        op.(name) = feval(op.(name), X, Y, Z);
    elseif ( isa(op.(name), 'surfacefun') )
        op.(name) = reshape([op.(name).vals{:}], [npts numPatches]);
    elseif ( isscalar(op.(name)) )
        op.(name) = repmat(op.(name), [1 numPatches]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT RHS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate non-constant RHSs if required:
if ( isa(rhs, 'function_handle') )
    rhs = feval(rhs, X(ii,:), Y(ii,:), Z(ii,:));
elseif ( isa(rhs, 'surfacefun') )
    vals = rhs.vals;
    rhs = reshape([vals{:}], [npts numPatches]);
    rhs = rhs(ii,:);
elseif ( isnumeric(rhs) && isscalar(rhs) )
    rhs = repmat(rhs, numIntPts, numPatches);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% SOLVE LOCAL PROBLEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
L = cell(numPatches, 1);

% Loop over each patch:
for k = 1:numPatches

    % Define the left and right edges for this patch:
    x = dom.x{k};
    y = dom.y{k};
    z = dom.z{k};
    edges = [ x(1) y(1) z(1) x(n)   y(n)   z(n)   nskel ;  % "Left" side
              x(1) y(1) z(1) x(end) y(end) z(end) nskel ;  % "Down" side
              x(n) y(n) z(n) x(end) y(end) z(end) nskel ]; % "Hypot" side

    A = zeros(npts);
    Dx = ux(:,k).*Du + vx(:,k).*Dv;
    Dy = uy(:,k).*Du + vy(:,k).*Dv;
    Dz = uz(:,k).*Du + vz(:,k).*Dv;
    J = dom.J{k}(:);

    if ( dom.singular(k) )

        % Assemble matrix:
        if ( flags.dxx ), A = A + op.dxx(:,k).*(J.*(Dx*Dx)-(Dx*J).*Dx); end
        if ( flags.dyy ), A = A + op.dyy(:,k).*(J.*(Dy*Dy)-(Dy*J).*Dy); end
        if ( flags.dzz ), A = A + op.dzz(:,k).*(J.*(Dz*Dz)-(Dz*J).*Dz); end
        if ( flags.dxy ), A = A + op.dxy(:,k).*(J.*(Dx*Dy)-(Dx*J).*Dy); end
        if ( flags.dyx ), A = A + op.dyx(:,k).*(J.*(Dy*Dx)-(Dy*J).*Dx); end
        if ( flags.dyz ), A = A + op.dyz(:,k).*(J.*(Dy*Dz)-(Dy*J).*Dz); end
        if ( flags.dzy ), A = A + op.dzy(:,k).*(J.*(Dz*Dy)-(Dz*J).*Dy); end
        if ( flags.dxz ), A = A + op.dxz(:,k).*(J.*(Dx*Dz)-(Dx*J).*Dz); end
        if ( flags.dzx ), A = A + op.dzx(:,k).*(J.*(Dz*Dx)-(Dz*J).*Dx); end
        if ( flags.dx  ), A = A + op.dx(:,k).*J.^2.*Dx;                 end
        if ( flags.dy  ), A = A + op.dy(:,k).*J.^2.*Dy;                 end
        if ( flags.dz  ), A = A + op.dz(:,k).*J.^2.*Dz;                 end
        if ( flags.b   ), A = A + op.b(:,k).*J.^3.*II;                  end

        % Construct solution operator:
        dA = decomposition(A(ii,ii), 'cod');
        Ainv = @(u) dA \ (J(ii).^3.*u);
        S = dA \ ([-A(ii,ee), J(ii).^3.*rhs(:,k)]);

        dx = L2S * (J(ee).^2.*Dx(ee,:));
        dy = L2S * (J(ee).^2.*Dy(ee,:));
        dz = L2S * (J(ee).^2.*Dz(ee,:));

        % The D2N map needs to be scaled on each side (e.g. when being
        % merged) to account for the Jacobian scaling which has been
        % factored out of the coordinate derivative maps. This scaling
        % is not known until the merge stage, as it depends on the
        % scaling of the neighboring patch.
        Jss = L2S * J(ee).^3;
        ItI_scl = {Jss(leftSkel); Jss(downSkel); Jss(hypotSkel)};

    else

        % Assemble matrix:
        if ( flags.dxx ), A = A + op.dxx(:,k).*(Dx*Dx); end
        if ( flags.dyy ), A = A + op.dyy(:,k).*(Dy*Dy); end
        if ( flags.dzz ), A = A + op.dzz(:,k).*(Dz*Dz); end
        if ( flags.dxy ), A = A + op.dxy(:,k).*(Dx*Dy); end
        if ( flags.dyx ), A = A + op.dyx(:,k).*(Dy*Dx); end
        if ( flags.dyz ), A = A + op.dyz(:,k).*(Dy*Dz); end
        if ( flags.dzy ), A = A + op.dzy(:,k).*(Dz*Dy); end
        if ( flags.dxz ), A = A + op.dxz(:,k).*(Dx*Dz); end
        if ( flags.dzx ), A = A + op.dzx(:,k).*(Dz*Dx); end
        if ( flags.dx  ), A = A + op.dx(:,k).*Dx;       end
        if ( flags.dy  ), A = A + op.dy(:,k).*Dy;       end
        if ( flags.dz  ), A = A + op.dz(:,k).*Dz;       end
        if ( flags.b   ), A = A + op.b(:,k).*II;        end

        dir_left  = II(ee_left,:);  neu_left  = NL(:,1,k).*Dx(ee_left,:)  + NL(:,2,k).*Dy(ee_left,:)  + NL(:,3,k).*Dz(ee_left,:);
        dir_down  = II(ee_down,:);  neu_down  = ND(:,1,k).*Dx(ee_down,:)  + ND(:,2,k).*Dy(ee_down,:)  + ND(:,3,k).*Dz(ee_down,:);
        dir_hypot = II(ee_hypot,:); neu_hypot = NH(:,1,k).*Dx(ee_hypot,:) + NH(:,2,k).*Dy(ee_hypot,:) + NH(:,3,k).*Dz(ee_hypot,:);

        bc_left  = neu_left  + 1i*eta*dir_left;
        bc_down  = neu_down  + 1i*eta*dir_down;
        bc_hypot = neu_hypot + 1i*eta*dir_hypot;

        leftLeaf  = 1:n;
        downLeaf  = [1 n+1:2:3*n-3];
        hypotLeaf = [n:2:3*n-4 3*n-3];
        corners = [1 n 3*n-3];

        bc = zeros(numBdyPts, npts);
        bc(leftLeaf,:)  = bc(leftLeaf,:)  + bc_left;
        bc(downLeaf,:)  = bc(downLeaf,:)  + bc_down;
        bc(hypotLeaf,:) = bc(hypotLeaf,:) + bc_hypot;
        bc(corners,:) = bc(corners,:)/2;

        BC = [ zeros(numIntPts, numBdyPts) ; eye(numBdyPts) ];
        RHS = [rhs(:,k); zeros(numBdyPts, 1)];

        LL = [ A(ipde,:) ; bc ];
        dA = matlab.internal.decomposition.DenseLU(LL);
        Ainv = @(u) solve(dA, u, false);
        S = Ainv([BC RHS]);

        ItI_scl = ItI_scl0;
    end

    % Append boundary points to solution operator and extract the
    % particular solution to store separately:
    u_part = S(:,end);
    S = S(:,1:end-1) * S2L;

    % Construct outgoing impedance operator:
    normal_d = [ neu_left  - 1i*eta*dir_left  ;
                 neu_down  - 1i*eta*dir_down  ;
                 neu_hypot - 1i*eta*dir_hypot ];
    normal_d = B2S * normal_d;

    % Construct the ItI map and particular flux:
    ItI = normal_d * S;
    du_part = normal_d * u_part;

    JJ = L2S * sqrt(J(ee));
    ww = wskel .* JJ;
    xyz = L2S * [x(ee) y(ee) z(ee)];

    % Assemble the patch:
    L{k} = surfaceop.leaf(dom, n, k, S, ItI, ItI_scl, u_part, du_part, edges, xyz, ww, Ainv, normal_d);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nl, nd, nh] = binormals(dom)
%BINORMALS   Compute the binormal vectors for a surfacemesh.

n = order(dom)+1;
eleft  = 1:n;
edown  = cumsum([1 (n:-1:2)]);
ehypot = cumsum([n (n-1:-1:1)]);

npts = size(dom.x{1}, 1);
sz = [npts 1 length(dom)];
xu = reshape([dom.xu{:}], sz); xv = reshape([dom.xv{:}], sz);
yu = reshape([dom.yu{:}], sz); yv = reshape([dom.yv{:}], sz);
zu = reshape([dom.zu{:}], sz); zv = reshape([dom.zv{:}], sz);

% Normal vectors to the surface (unnormalized)
nl = -[xu(eleft,:,:)  yu(eleft,:,:)  zu(eleft,:,:)];
nd = -[xv(edown,:,:)  yv(edown,:,:)  zv(edown,:,:)];
nh =  [xu(ehypot,:,:) yu(ehypot,:,:) zu(ehypot,:,:)] + ...
      [xv(ehypot,:,:) yv(ehypot,:,:) zv(ehypot,:,:)];

% Tangent vectors to the element boundary (normalized)
tl = normalize([xv(eleft,:,:)   yv(eleft,:,:)   zv(eleft,:,:)]);
td = normalize([xu(edown,:,:)   yu(edown,:,:)   zu(edown,:,:)]);
th = normalize( -[xv(ehypot,:,:) yv(ehypot,:,:) zv(ehypot,:,:)] + ...
                 [xu(ehypot,:,:) yu(ehypot,:,:) zu(ehypot,:,:)] );

% Binormal vectors (normalized)
nl = normalize(nl - tl .* sum(nl.*tl, 2));
nd = normalize(nd - td .* sum(nd.*td, 2));
nh = normalize(nh - th .* sum(nh.*th, 2));

end

function v = normalize(v)

v = v ./ sqrt(v(:,1,:).^2 + v(:,2,:).^2 + v(:,3,:).^2);

end

function P = skel2leaf(nleaf, nskel)
%SKEL2LEAF   Boundary interpolation matrix.
%   SKEL2LEAF(NLEAF, NSKEL) returns the (3*NLEAF-3) x 3*NSKEL matrix that
%   maps 3 pieces of length-NSKEL first-kind boundary values to 3*NLEAF-3
%   second-kind boundary values, including the corners. At each corner, the
%   average of the two interpolated values is used.

[xskel, ~, wskel] = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(nleaf, 2);
B = barymat(xleaf, xskel, wskel);

% Skeleton indices for each side
leftSkel  = 1:nskel;
downSkel  = nskel+1:2*nskel;
hypotSkel = 2*nskel+1:3*nskel;

% Leaf indices for each side
leftLeaf  = 1:nleaf;
downLeaf  = [1 nleaf+1:2:3*nleaf-3];
hypotLeaf = [nleaf:2:3*nleaf-4 3*nleaf-3];

P = zeros(3*nleaf-3, 3*nskel);
P(leftLeaf,  leftSkel)  = B;
P(downLeaf,  downSkel)  = B;
P(hypotLeaf, hypotSkel) = B;

% Average the corners:
corners = [1 nleaf 3*nleaf-3];
P(corners,:) = P(corners,:)/2;

end

function P = leaf2skel(nskel, nleaf)
%LEAF2SKEL   Boundary interpolation matrix.
%   LEAF2SKEL(NSKEL, NLEAF) returns the 3*NSKEL x (3*NLEAF-3) matrix that
%   maps 3*NLEAF-3 second-kind boundary values to 3 pieces of length-NSKEL
%   first-kind boundary values.

[xskel, ~, wskel] = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(nleaf, 2);
B = barymat(xskel, xleaf, wleaf);

% Skeleton indices for each side
leftSkel  = 1:nskel;
downSkel  = nskel+1:2*nskel;
hypotSkel = 2*nskel+1:3*nskel;

% Leaf indices for each side
leftLeaf  = 1:nleaf;
downLeaf  = [1 nleaf+1:2:3*nleaf-3];
hypotLeaf = [nleaf:2:3*nleaf-4 3*nleaf-3];

P = zeros(3*nskel, 3*nleaf-3);
P(leftSkel,  leftLeaf)  = B;
P(downSkel,  downLeaf)  = B;
P(hypotSkel, hypotLeaf) = B;

end

function P = bdy2skel(nskel, nleaf)
%BDY2SKEL   Boundary interpolation matrix.
%   BDY2SKEL(NSKEL, NLEAF) returns the 3*NSKEL x 3*NLEAF matrix that maps
%   3 pieces of length-NLEAF second-kind boundary values to 3 pieces of
%   length-NSKEL first-kind boundary values.

[xskel, ~, wskel] = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(nleaf, 2);
B = barymat(xskel, xleaf, wleaf);

% Skeleton indices for each side
leftSkel  = 1:nskel;
downSkel  = nskel+1:2*nskel;
hypotSkel = 2*nskel+1:3*nskel;

% Leaf boundary indices for each side
leftBdy  = 1:nleaf;
downBdy  = nleaf+1:2*nleaf;
hypotBdy = 2*nleaf+1:3*nleaf;

P = zeros(3*nskel, 3*nleaf);
P(leftSkel,  leftBdy)  = B;
P(downSkel,  downBdy)  = B;
P(hypotSkel, hypotBdy) = B;

end
