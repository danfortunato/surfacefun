function [i1, i2, i4a, i4b, dom, edgesAB] = intersect(a, b)
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
tol = 1e-12 * scl;
% Must check in both directions:
[iA1, iB1] = intersectTol(edgesA(:,[1 2 3 4 5 6]), edgesB(:,1:6), tol);
[iA2, iB2] = intersectTol(edgesA(:,[4 5 6 1 2 3]), edgesB(:,1:6), tol);
iA = [iA1 ; iA2];
iB = [iB1 ; iB2];
assert(numel(iA) == numel(iB), 'Intersection failed.');

% Determine indices corresponding to intersecting edges:
pA = edgesA(:,7) - 2;
ppA = cumsum([0 ; pA]);
i4a = [];
for k = 1:numel(iA)
    i4a = [i4a ; ppA(iA(k)) + (1:pA(iA(k))).'];
end
pB = edgesB(:,7) - 2;
ppB = cumsum([0 ; pB]);
i4b = [];
for k = 1:numel(iB)
    i4b = [i4b ; ppB(iB(k)) + (1:pB(iB(k))).'];
end

% i1 and i2 are remaining points (i.e., those not in the intersection).
i1 = (1:sum(pA)).'; i1(i4a) = [];
i2 = (1:sum(pB)).'; i2(i4b) = [];

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
