function [i1, i2, i4a, i4b, flip1, flip2, scl1, scl2, sclAB, dom, edgesAB] = intersect(a, b)
%INTERSECT   Compute the indices of the glue between two patches.
%   [I1, I2, I4A, I4B, L2G1, L2G2, SCL1, SCL2, SCLAB] = INTERSECT(A, B)
%   returns the indices of the glue w.r.t. A.edges and B.edges of two
%   patches A and B. For consistency with the paper by Martinsson:
%
%       I1 : Indices of DOFs on A.edges which are not on B.edges
%       I2 : Indices of DOFs on B.edges which are not on A.edges
%       I4A: Indices of DOFs on A.edges which are in the intersection
%       I4B: Indices of DOFs on B.edges which are in the intersection
%
%   If the boundaries being merged contain flipped regions (so that DOFs
%   are ordered differently local to each patch) then FLIP1 and FLIP2
%   encode how to map from local DOFs in A and B to global DOFs along
%   the shared interface. The matrices SCL1 and SCL2 are scale vectors for
%   the Jacobians that have been factored out of the Dirichlet-to-Neumann
%   maps for patches A and B. SCLAB is a cell array of scalars and/or
%   function handles defining those algebraic expressions for the parent's
%   edges.
%
%   [I1, I2, I4A, I4B, DOM, EDGESAB] = INTERSECT(A, B) returns also the
%   domain and the edges of the intersection.

% The new domain will be NaN except in the special case when we are
% merging two rectangluar patches horizontally or vertically.
dom = NaN;

% Above to eventually be replaced by:
edgesA = a.edges;
edgesB = b.edges;

sclx = 0;
scly = 0;
sclz = 0;
for k = 1:4
    sclx = max(sclx, abs(edgesA(k,4) - edgesA(k,1)));
    sclx = max(sclx, abs(edgesB(k,4) - edgesB(k,1)));
    scly = max(scly, abs(edgesA(k,5) - edgesA(k,2)));
    scly = max(scly, abs(edgesB(k,5) - edgesB(k,2)));
    sclz = max(sclz, abs(edgesA(k,6) - edgesA(k,3)));
    sclz = max(sclz, abs(edgesB(k,6) - edgesB(k,3)));
end
scl = max([sclx scly sclz]);

% Check for intersecting edges (with a tolerance):
tol = 1e-8 * scl;
% Must check in both directions:
digits = ceil(abs(log10(tol)));
%[~, iA1, iB1] = intersect(round(edgesA(:,[1 2 3 4 5 6]), digits), ...
%                          round(edgesB(:,1:6), digits), 'rows');
[iA1, iB1] = intersectTol(edgesA(:,[1 2 3 4 5 6]), edgesB(:,1:6), tol);
[iA2, iB2] = intersectTol(edgesA(:,[4 5 6 1 2 3]), edgesB(:,1:6), tol);
iA = [iA1 ; iA2];
iB = [iB1 ; iB2];
assert(numel(iA) == numel(iB), 'Intersection failed.');

% Determine indices corresponding to intersecting edges:
pA = edgesA(:,7);
ppA = cumsum([0 ; pA]);
i4a = [];
for k = 1:numel(iA)
    i4a = [i4a ; ppA(iA(k)) + (1:pA(iA(k))).'];
end
pB = edgesB(:,7);
ppB = cumsum([0 ; pB]);
i4b = [];
for k = 1:numel(iB)
    i4b = [i4b ; ppB(iB(k)) + (1:pB(iB(k))).'];
end

flip1 = speye(length(i4a));
flip2 = cell(numel(iB), 1);
for k = 1:numel(iB1)
    flip2{k} = speye(pB(iB1(k)));
end
for k = 1:numel(iB2)
    flip2{numel(iB1)+k} = fliplr(speye(pB(iB2(k))));
end
if ( ~isempty(flip2) )
    flip2 = matlab.internal.math.blkdiag(flip2{:});
else
    flip2 = flip1;
end

% i1 and i2 are remaining points (i.e., those not in the intersection).
i1 = (1:sum(pA)).'; i1(i4a) = [];
i2 = (1:sum(pB)).'; i2(i4b) = [];

% TODO: Deal with flipping?
%flip1 = ones(size(i4a)); flip1 = diag(flip1);
%flip2 = ones(size(i4b)); flip2 = diag(flip2);
% flip1 = ones(size(i4a));
% flip2 = ones(sum(pB(iB1)),1);
% for k = 1:numel(iB2)
%     e = ones(pB(iB2(k)),1); e(2:2:end) = -1;
%     flip2 = [flip2 ; e];
% end
% if ( isempty(flip1) )
%     % Force 0x1 rather than empty. (Not sure why this is required.)
%     flip1 = ones(0,1);
%     flip2 = ones(0,1);
% end

%% Construct operators for p-adaptivity and Jacobian scaling.

sclA = a.D2N_scl;
sclB = b.D2N_scl;
scl1 = cat(1, sclA{iA});
scl2 = cat(1, sclB{iB});
if ( isempty(scl1) )
    scl1 = zeros(0,1);
    scl2 = zeros(0,1);
end

% Concatenate the scaling functions for the parent's edges.
if ( ~isempty(iA) )
    sclA(iA) = [];
    sclB(iB) = [];
end
sclAB = [sclA ; sclB];

edgesA(iA,:) = [];
edgesB(iB,:) = [];
edgesAB = [edgesA ; edgesB];

end

function [AI, BI] = intersectTol(A, B, tol)
% Borrowed from https://www.mathworks.com/matlabcentral/...
% answers/444501-how-to-use-intersect-command-with-a-tolerance-value
   n = size(A,1);
   M  = zeros(n,1);
   % Collect the index of the first occurrence in B for every A:
   for k = 1:n
      dist = sum(abs(A(k,:) - B), 2);    % 1-norm
      idx  = find(dist < tol, 1);        % Absolute tolerance
      % iddx = find(dist ./ A(iA) < tol, 1);  % Relative tolerance
      if ( ~isempty(idx) )
         M(k) = idx;
      end
   end
   AI = find(M);
   BI = M(AI); 
end
