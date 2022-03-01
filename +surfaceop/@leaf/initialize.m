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
ii([1,end],[1,end]) = true;       % Treat corners as interior.
ee = ~ii;                         % Boundary indices.
leftIdx  = 1:n-2;
rightIdx = n-1:2*(n-2);
downIdx  = 2*(n-2)+1:3*(n-2);
upIdx    = 3*(n-2)+1:4*(n-2);
numBdyPts = sum(ee(:)); 
numIntPts = sum(ii(:));
ibc = 3*(n-2)+1;
ss = [1:n-2, ibc:4*(n-2), n-1:2:ibc-1, n:2:ibc-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = diffmat(n);
I = eye(n);
II = kron(I, I);
Du = kron(D, I);
Dv = kron(I, D);

% Interpolation operator for corner values:
Xii = X(1,2:(n-1)).';
B = [Xii-1, -1-Xii].'; B(:,1:2:end) = -B(:,1:2:end);
if ( mod(n-1, 2) )
    B(2,:) = -B(2,:);
end

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
    edges = [ x(1,1) y(1,1) z(1,1) x(n,1) y(n,1) z(n,1) n ;  % "Left" side
              x(1,n) y(1,n) z(1,n) x(n,n) y(n,n) z(n,n) n ;  % "Right" side
              x(1,1) y(1,1) z(1,1) x(1,n) y(1,n) z(1,n) n ;  % "Down" side
              x(n,1) y(n,1) z(n,1) x(n,n) y(n,n) z(n,n) n ]; % "Up" side

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

    %[Dx, Dy, Dz] = diffs(x, y, z);
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

    for name = fieldnames(op).'
        name = name{1};
        if ( isa(op.(name), 'function_handle') )
            op.(name) = feval(op.(name), x, y, z);
        elseif ( isa(op.(name), 'surfacefun') )
            op.(name) = op.(name).vals{k};
        end
    end

    A = zeros(n^2);
    if ( op.dxx ~= 0 ), A = A + op.dxx(:).*Dxx; end
    if ( op.dyy ~= 0 ), A = A + op.dyy(:).*Dyy; end
    if ( op.dzz ~= 0 ), A = A + op.dzz(:).*Dzz; end
    if ( op.dxy ~= 0 ), A = A + op.dxy(:).*Dxy; end
    if ( op.dyz ~= 0 ), A = A + op.dyz(:).*Dyz; end
    if ( op.dxz ~= 0 ), A = A + op.dxz(:).*Dxz; end
    if ( op.dx  ~= 0 ), A = A + op.dx(:).*Dx;   end
    if ( op.dy  ~= 0 ), A = A + op.dy(:).*Dy;   end
    if ( op.dz  ~= 0 ), A = A + op.dz(:).*Dz;   end
    if ( op.b   ~= 0 ), A = A + op.b(:).*II;    end

    % Construct solution operator:
    %Ainv = @(u) A(ii,ii) \ u;
    if ( dom.singular(k) )
        dA = decomposition(A(ii,ii), 'cod');
        %dA = decomposition(A(ii,ii));
        Ainv = @(u) dA \ (dom.J{k}(ii).^3 .* u);
        %S = Ainv([-A(ii,ee), rhs_eval]);
        S = dA \ ([-A(ii,ee), dom.J{k}(ii).^3.*rhs_eval]);
    else
        dA = decomposition(A(ii,ii));
        Ainv = @(u) dA \ u;
        %A1 = inv(A(ii,ii));
        %Ainv = @(u) A1 * u;
        %[U1, S1, V1] = svd(A(ii,ii));
        %U1 = U1'./diag(S1);
        %Ainv = @(u) V1*(U1*u);
        S = Ainv([-A(ii,ee), rhs_eval]);
    end

    % Replace solution operator for corners with interp conditions:
    S([1:2,end-1:end],:) = 0;
    S(1:2,1:n-2) = B;
    S([end-1,end],end-n+2:end-1) = B;

    % Append boundary points to solution operator:
    tmpS = zeros(n^2, size(S, 2));
    tmpS(ii,:) = S;
    tmpS(ee,:) = eye(numBdyPts, numBdyPts+1);
    S = tmpS;
    S = S(:,[ss end]);

    % Construct the D2N map:
    dx = Dx(ee,:) * S; dx = dx(ss,:);
    dy = Dy(ee,:) * S; dy = dy(ss,:);
    dz = Dz(ee,:) * S; dz = dz(ss,:);
    [nl, nr, nd, nu] = normals(x, y, z);
    D2N = zeros(numBdyPts, numBdyPts+1);
    D2N(leftIdx,:)  = nl(:,1).*dx(leftIdx,:)  + nl(:,2).*dy(leftIdx,:)  + nl(:,3).*dz(leftIdx,:);
    D2N(rightIdx,:) = nr(:,1).*dx(rightIdx,:) + nr(:,2).*dy(rightIdx,:) + nr(:,3).*dz(rightIdx,:);
    D2N(downIdx,:)  = nd(:,1).*dx(downIdx,:)  + nd(:,2).*dy(downIdx,:)  + nd(:,3).*dz(downIdx,:);
    D2N(upIdx,:)    = nu(:,1).*dx(upIdx,:)    + nu(:,2).*dy(upIdx,:)    + nu(:,3).*dz(upIdx,:);

    dx = Dx(ee,:); dx = dx(ss,:);
    dy = Dy(ee,:); dy = dy(ss,:);
    dz = Dz(ee,:); dz = dz(ss,:);
    normal_d = zeros(numBdyPts, n^2);
    normal_d(leftIdx,:)  = nl(:,1).*dx(leftIdx,:)  + nl(:,2).*dy(leftIdx,:)  + nl(:,3).*dz(leftIdx,:);
    normal_d(rightIdx,:) = nr(:,1).*dx(rightIdx,:) + nr(:,2).*dy(rightIdx,:) + nr(:,3).*dz(rightIdx,:);
    normal_d(downIdx,:)  = nd(:,1).*dx(downIdx,:)  + nd(:,2).*dy(downIdx,:)  + nd(:,3).*dz(downIdx,:);
    normal_d(upIdx,:)    = nu(:,1).*dx(upIdx,:)    + nu(:,2).*dy(upIdx,:)    + nu(:,3).*dz(upIdx,:);

    % The D2N map needs to be scaled on each side (e.g. when being
    % merged) to account for the Jacobian scaling which has been
    % factored out of the coordinate derivative maps. This scaling
    % is not known until the merge stage, as it depends on the
    % scaling of the neighboring patch.
    if ( dom.singular(k) )
        D2N_scl = cell(4, 1);
        J = dom.J{k}.^3; Jee = J(ee); Jss = Jee(ss);
        D2N_scl{1} = Jss(leftIdx); % Left
        D2N_scl{2} = Jss(rightIdx);  % Right
        D2N_scl{3} = Jss(downIdx); % Down
        D2N_scl{4} = Jss(upIdx);  % Up
    else
        D2N_scl = {ones(n-2,1), ones(n-2,1), ones(n-2,1), ones(n-2,1)}.';
    end

    % Extract the particular solution to store separately:
    u_part = S(:,end); S = S(:,1:end-1);
    du_part = D2N(:,end); D2N = D2N(:,1:end-1);

    % Assemble the patch:
    xee = x(ee);
    yee = y(ee);
    zee = z(ee);
    xyz = [xee(ss) yee(ss) zee(ss)];
    L{k} = surfaceop.leaf(dom, k, S, D2N, D2N_scl, u_part, du_part, edges, xyz, Ainv, normal_d);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Dx, Dy, Dz] = diffs(x, y, z)
%DIFFS   Differential operators on surfaces.
%   [DX, DY, DZ] = DIFFS(X, Y, Z) returns the tangential derivative
%   operators DX, DY, and DZ over the surface defined by the coordinates
%   (X, Y, Z). The coordinates must be given as N x N matrices sampled at
%   tensor-product Chebyshev nodes. DX, DY, and DZ are operators of size
%   N^2 x N^2.

persistent D
n = size(x, 1);
if ( size(D, 1) ~= n )
    D = diffmat(n);
end

xu = x * D.'; xv = D * x;
yu = y * D.'; yv = D * y;
zu = z * D.'; zv = D * z;
E = xu.*xu + yu.*yu + zu.*zu;
G = xv.*xv + yv.*yv + zv.*zv;
F = xu.*xv + yu.*yv + zu.*zv;
J = E.*G - F.^2;
ux = (G.*xu-F.*xv)./J; vx = (E.*xv-F.*xu)./J;
uy = (G.*yu-F.*yv)./J; vy = (E.*yv-F.*yu)./J;
uz = (G.*zu-F.*zv)./J; vz = (E.*zv-F.*zu)./J;
I = eye(n); % or speye?
Du = kron(D, I);
Dv = kron(I, D);
Dx = ux(:).*Du + vx(:).*Dv;
Dy = uy(:).*Du + vy(:).*Dv;
Dz = uz(:).*Du + vz(:).*Dv;

end

function [nl, nr, nd, nu] = normals(x, y, z)
%NORMALS   Outward pointing normal vectors to the edges of a mapping.

persistent D
n = size(x, 1);
if ( size(D, 1) ~= n )
    D = diffmat(n);
end

xu = x * D.'; xv = D * x;
yu = y * D.'; yv = D * y;
zu = z * D.'; zv = D * z;

nl = -normalize([xu(:,1)   yu(:,1)   zu(:,1)]);
nr =  normalize([xu(:,n)   yu(:,n)   zu(:,n)]);
nd = -normalize([xv(1,:).' yv(1,:).' zv(1,:).']);
nu =  normalize([xv(n,:).' yv(n,:).' zv(n,:).']);

tangent = normalize([xv(:,1) yv(:,1) zv(:,1)]);
nl = nl - tangent .* dot(nl, tangent, 2);
tangent = normalize([xv(:,n) yv(:,n) zv(:,n)]);
nr = nr - tangent .* dot(nr, tangent, 2);
tangent = normalize([xu(1,:).' yu(1,:).' zu(1,:).']);
nd = nd - tangent .* dot(nd, tangent, 2);
tangent = normalize([xu(n,:).' yu(n,:).' zu(n,:).']);
nu = nu - tangent .* dot(nu, tangent, 2);

nl = normalize(nl);
nr = normalize(nr);
nd = normalize(nd);
nu = normalize(nu);

nl([1,n],:) = [];
nr([1,n],:) = [];
nd([1,n],:) = [];
nu([1,n],:) = [];

end

function v = normalize(v)

v = v ./ sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);

end
