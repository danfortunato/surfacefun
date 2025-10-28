function mesh(dom, varargin)
%MESH   Plot the mesh of a surface.

surfaceStyle = {'FaceColor',        'w',     ...
                'AmbientStrength',   0.6,    ...
                'DiffuseStrength',   0.4,    ...
                'SpecularStrength',  0.3 };

holdState = ishold();

x = dom.x;
y = dom.y;
z = dom.z;
vn = dom.facenormals;
ne = length(dom);

if ( ~all(dom.ptype == dom.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

switch ( dom.ptype(1) )
    case 'tri'
        n = order(dom)+1;
        npts = length(dom.x{1});
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

        hold on
        patch('Faces', T_all, ...
              'Vertices', [x_all y_all z_all], ...
              'FaceVertexCData', 0*x_all, ...
              surfaceStyle{:}, ...
              varargin{:});

    case 'quad'
        for k = 1:length(dom)
            surf(x{k}, y{k}, z{k}, 0*x{k}, 'FaceColor', 'None', varargin{:})
            hold on
            % If the plot is not being added to another then plot the surface so
            % that the lines are more easily discernable.
            if ( ~holdState )
                % Plot the surface, making it slightly smaller so lines show up
                % more clearly.
                scl = 0.01;
                surface(x{k} - scl*vn{k}(:,:,1), ...
                        y{k} - scl*vn{k}(:,:,2), ...
                        z{k} - scl*vn{k}(:,:,3), ...
                        0*x{k}, ...
                        'EdgeColor', 'none', ...
                        surfaceStyle{:});
            end
        end
end

if ( ~holdState )
    view(3)
    axis equal
    hold off
end

end
