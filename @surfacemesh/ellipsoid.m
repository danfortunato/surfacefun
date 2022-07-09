function dom = ellipsoid(n, abc, varargin)
%ELLIPSOID   Create a cubed ellipsoid mesh.

if ( nargin < 2 )
    abc = [1 1 1];
end

% Make a sphere:
dom = surfacemesh.sphere(n, varargin{:});

% Stretch it by (a, b, c):
[x, y, z] = deal(dom.x, dom.y, dom.z);
for k = 1:length(dom)
    x{k} = abc(1) * x{k};
    y{k} = abc(2) * y{k};
    z{k} = abc(3) * z{k};
end
dom = surfacemesh(x, y, z);

end
