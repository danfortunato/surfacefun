function varargout = mesh(f, varargin)
%MESH   3-D mesh surface of a SURFACEFUN.
%   MESH(F) plots the tensor-product Chebyshev grids for the patches of F,
%   colored according to the values of F.
%
%   See also PLOT, SURF, CONTOUR.

surfaceStyle = {'FaceColor',        'w',     ...
                'AmbientStrength',   0.6,    ...
                'DiffuseStrength',   0.4,    ...
                'SpecularStrength',  0.3 };

holdState = ishold();

if ( ~all(f.domain.ptype == f.domain.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

switch ( f.domain.ptype(1) )
    case 'tri'
        x = f.domain.x;
        y = f.domain.y;
        z = f.domain.z;
        ne = length(f.domain);
        n = order(f.domain)+1;
        npts = length(f.domain.x{1});
        T = trilattice(n);
        ntri = size(T, 1);
        T_all = repmat(T, [1 1 ne]) + reshape((0:ne-1)*npts, 1, 1, []);
        T_all = permute(T_all, [1 3 2]);
        T_all = reshape(T_all, ntri*ne, 3);
        x = [x{:}];
        y = [y{:}];
        z = [z{:}];
        x_all = x(:);
        y_all = y(:);
        z_all = z(:);
        f_all = reshape([f.vals{:}], [], 1);

        hold on
        patch('Faces', T_all, ...
              'Vertices', [x_all y_all z_all], ...
              'FaceVertexCData', f_all, ...
              'EdgeColor', 'interp', ...
              surfaceStyle{:}, ...
              varargin{:});

    case 'quad'
        [x, y, z] = getGrid(f);
        for k = 1:length(f)
            u = f.vals{k};
            if ( ~isreal(u) )
                u = abs(u);
            end
            h(k) = mesh(x{k}, y{k}, z{k}, u, varargin{:}); %#ok<AGROW>
            hold on
        end
end

if ( ~holdState )
    view(3)
    axis equal
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
