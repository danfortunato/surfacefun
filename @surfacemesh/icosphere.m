function dom = icosphere(n, nref)
%ICOSPHERE   Create an icosahedral sphere mesh.

if ( nargin < 2 )
    nref = 0;
end

% Create a regular unit icosahedron
[v, f] = icosahedron();

% Recursively subdivide triangle faces
for gen = 1:nref
    nf = size(f, 1);
    f_ = zeros(4*nf, 3);
    % Loop over each triangle
    for i = 1:nf
        tri = f(i,:);
        % Calculate mid points (add new points to v)
        [a, v] = getMidPoint(tri(1), tri(2), v);
        [b, v] = getMidPoint(tri(2), tri(3), v);
        [c, v] = getMidPoint(tri(3), tri(1), v);
        % Generate new subdivision triangles
        nfc = [ tri(1) a c ;
                tri(2) b a ;
                tri(3) c b ;
                a      b c ];
        % Replace triangle with subdivision
        idx = 4*(i-1)+1:4*i;
        f_(idx,:) = nfc;
    end
    f = f_;
end

% Remove duplicate vertices
[v, ~, ix] = unique(v, 'rows');

% Reassign faces to trimmed vertex list and remove any duplicate faces
f = unique(ix(f), 'rows');

dom0 = [0 0; 0 1; 1 0];
dom0 = [dom0 zeros(3,1)];
[uu0, vv0] = trianglepts(n, type='cheb2');
uvw0 = [uu0 vv0 zeros(size(uu0)) ones(size(uu0))].';

nf = size(f, 1);
x = cell(nf, 1);
y = cell(nf, 1);
z = cell(nf, 1);
for k = 1:nf
    tri = f(k,:);
    dom = v(tri,:);
    A = affine(dom0, dom);
    xyz = A*uvw0;
    nrm = sqrt(xyz(1,:).^2 + xyz(2,:).^2 + xyz(3,:).^2);
    x{k} = xyz(1,:).' ./ nrm.';
    y{k} = xyz(2,:).' ./ nrm.';
    z{k} = xyz(3,:).' ./ nrm.';
end

dom = surfacemesh(x, y, z, 'tri');

end

function A = affine(from, to)
%AFFINE   Affine transformation.
%   A = AFFINE(FROM, TO) 

if ( ~ismatrix(from) || ~ismatrix(to) || ~all(size(from) == size(to)) )
    error('Invalid points.');
end

n   = size(from, 1);
dim = size(from, 2);

w = ones(1, n);
A = [to.'; w] / [from.' ; w];

end

function [i, v] = getMidPoint(t1, t2, v)
%GETMIDPOINT   Calculate midpoint between two vertices.
%   Calculate new vertex in sub-division and normalise to unit length
%   then find or add it to v and return index.

% Get vertex positions
p1 = v(t1,:);
p2 = v(t2,:);

% Calculate mid point (on unit sphere)
pm = (p1 + p2) / 2;
pm = pm / norm(pm);

% Add to vertices list, return index
i = size(v, 1) + 1;
v = [v ; pm];

end

function [v,f] = icosahedron()
%ICOSAHEDRON   Create a unit regular icosahedron.
%   Returns 12 vertex and 20 face values.

t = (1 + sqrt(5)) / 2;

v = [ -1  t  0 ;
       1  t  0 ;
      -1 -t  0 ;
       1 -t  0 ;
       0 -1  t ;
       0  1  t ;
       0 -1 -t ;
       0  1 -t ;
       t  0 -1 ;
       t  0  1 ;
      -t  0 -1 ;
      -t  0  1 ];

v = v / norm(v(1,:));

f = [ 1 12  6 ;
      1  6  2 ;
      1  2  8 ;
      1  8 11 ;
      1 11 12 ;
      2  6 10 ;
      6 12  5 ;
     12 11  3 ;
     11  8  7 ;
      8  2  9 ;
      4 10  5 ;
      4  5  3 ;
      4  3  7 ;
      4  7  9 ;
      4  9 10 ;
      5 10  6 ;
      3  5 12 ;
      7  3 11 ;
      9  7  8 ;
     10  9  2 ];

end
