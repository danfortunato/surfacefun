function varargout = quiver(f, scl, m, varargin)
%QUIVER   Quiver plot of a SURFACEFUNV.
%   QUIVER(F) plots the vector field of the SURFACEFUNV F on its surface.
%
%   QUIVER(F, SCL) scales the vectors by the scalar SCL.
%
%   QUIVER(F, SCL, M) plots M x M vectors on each patch of the surface.
%
%   H = QUIVER(F) returns a handle to the figure window.

% Empty check:
if ( isempty(f) )
    h = quiver([]);
    if ( nargout > 0 )
        varargout = {h};
    end
    return
end

holdState = ishold();

if ( nargin < 2 )
    scl = 1;
end

if ( nargin < 3 )
    m = 10;
end

dom = f.components{1}.domain;
x = dom.x;
y = dom.y;
z = dom.z;
vn = dom.facenormals;

xx = zeros(m, m, length(dom)); fx = zeros(m, m, length(dom));
yy = zeros(m, m, length(dom)); fy = zeros(m, m, length(dom));
zz = zeros(m, m, length(dom)); fz = zeros(m, m, length(dom));
n = size(x{1}, 1);
B = bary(linspace(-1, 1, m).', eye(n));
for k = 1:length(dom)
    xx(:,:,k) = B * x{k} * B.';
    yy(:,:,k) = B * y{k} * B.';
    zz(:,:,k) = B * z{k} * B.';
    fx(:,:,k) = B * f.components{1}.vals{k} * B.';
    fy(:,:,k) = B * f.components{2}.vals{k} * B.';
    fz(:,:,k) = B * f.components{3}.vals{k} * B.';
    if ( ~holdState )
        r = 0.01;
        surface(x{k} - r*vn{k}(:,:,1), ...
                y{k} - r*vn{k}(:,:,2), ...
                z{k} - r*vn{k}(:,:,3), ...
                0*x{k}, 'FaceColor', 'w', 'EdgeColor', 'None');
        hold on
    end
end

h = quiver3(xx, yy, zz, scl*fx, scl*fy, scl*fz, 0, varargin{:});

if ( ~holdState )
    view(3)
    axis equal
    grid on
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
