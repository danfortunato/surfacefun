function L = initialize(dom, rhs)
%INITIALIZE   Initialize an array of LEAF objects.
%   L = SURFACESEM.LEAF.INITIALIZE(DOM) returns a cell array L of LEAF
%   objects which contain the solution and D2N operators for Poisson's
%   equation on the domain DOM with zero righthand side.
%
%   L = SURFACESEM.LEAF.INITIALIZE(DOM, RHS) is as above, but with the
%   righthand side RHS, which may be a scalar or a function handle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(isstruct(dom), 'Invalid domain.');

if ( nargin < 2 )
    % Default to homogeneous problem:
    rhs = 0;
end

numPatches = length(dom);
n = size(dom(1).x, 1);

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
    domk = dom(k);
    x = domk.x;
    y = domk.y;
    z = domk.z;
    edges = [ x(n,1) y(n,1) z(n,1) x(1,1) y(1,1) z(1,1) n ;  % "Left" side
              x(n,n) y(n,n) z(n,n) x(1,n) y(1,n) z(1,n) n ;  % "Right" side
              x(1,1) y(1,1) z(1,1) x(1,n) y(1,n) z(1,n) n ;  % "Down" side
              x(n,1) y(n,1) z(n,1) x(n,n) y(n,n) z(n,n) n ]; % "Up" side

    % Evaluate non-constant RHSs if required:
    if ( isa(rhs, 'function_handle') || isa(rhs, 'chebfun2') )
        rhs_eval = feval(rhs, x(ii), y(ii), z(ii));
    elseif ( iscell(rhs) )
        rhs_eval = rhs{k};
    else
        rhs_eval = rhs;
    end

    [Dx, Dy, Dz] = diffs(x, y, z);
    A = Dx^2 + Dy^2 + Dz^2;

    % Construct solution operator:
    S = A(ii,ii) \ [-A(ii,ee), rhs_eval];
    Ainv = @(u) A(ii,ii) \ u;

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

    % Assemble the patch:
    xee = x(ee);
    yee = y(ee);
    zee = z(ee);
    xyz = [xee(ss) yee(ss) zee(ss)];
    L{k} = SurfaceSEM.Leaf(domk, S, D2N, edges, xyz, Ainv);

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
