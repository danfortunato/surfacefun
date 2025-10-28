function wireframe(dom, varargin)
%WIREFRAME   Plot the wireframe of a surface.

parser = inputParser;
parser.KeepUnmatched = true;
parser.addParameter('surface', 'auto', @(s) contains(lower(s), {'auto', 'on', 'off'}));
parser.addParameter('edges',   'auto', @(s) contains(lower(s), {'auto', 'on', 'off'}));
parse(parser, varargin{:});
showSurface = parser.Results.surface;
showEdges   = parser.Results.edges;
varargin = namedargs2cell(parser.Unmatched);

defaultStyle = {'Color',     'k', ...
                'LineStyle', '-', ...
                'LineWidth',  1};

surfaceStyle = {'FaceColor',        'w',     ...
                'EdgeColor',        'none',  ...
                'AmbientStrength',   0.6,    ...
                'DiffuseStrength',   0.4,    ...
                'SpecularStrength',  0.3 };

holdState = ishold();

vn = dom.facenormals;
ne = length(dom);
scl = 0.005;

if ( all(dom.ptype == 'tri') )

    npts = length(dom.x{1});
    n = (sqrt(8*npts+1)-1) / 2;
    x = [dom.x{:}];
    y = [dom.y{:}];
    z = [dom.z{:}];
    e1 = 1:n;
    e2 = cumsum([1 (n:-1:2)]);
    e3 = cumsum([n (n-1:-1:1)]);
    X = [ x(e1,:) ; nan(1,ne) ; x(e2,:) ; nan(1,ne) ; x(e3,:) ; nan(1,ne) ];
    Y = [ y(e1,:) ; nan(1,ne) ; y(e2,:) ; nan(1,ne) ; y(e3,:) ; nan(1,ne) ];
    Z = [ z(e1,:) ; nan(1,ne) ; z(e2,:) ; nan(1,ne) ; z(e3,:) ; nan(1,ne) ];

    if ( ~strcmpi(showEdges, 'off') )
        plot3(X(:), Y(:), Z(:), defaultStyle{:}, varargin{:})
    end

    % If the plot is not being added to another then plot the surface so
    % that the lines are more easily discernable.
    if ( (~holdState && strcmpi(showSurface, 'auto')) || strcmpi(showSurface, 'on') )
        % Plot the surface, making it slightly smaller so lines show up
        % more clearly.

        T = trilattice(n);
        ntri = size(T, 1);
        T_all = repmat(T, [1 1 ne]) + reshape((0:ne-1)*npts, 1, 1, []);
        T_all = permute(T_all, [1 3 2]);
        T_all = reshape(T_all, ntri*ne, 3);
        vn_all = reshape([vn{:}], npts*ne, 3);
        x_all = x(:) - scl*vn_all(:,1);
        y_all = y(:) - scl*vn_all(:,2);
        z_all = z(:) - scl*vn_all(:,3);

        hold on
        patch('Faces',    T_all,               ...
              'Vertices', [x_all y_all z_all], ...
              'FaceVertexCData', 0*x_all,      ...
              surfaceStyle{:});
    end

elseif ( all(dom.ptype == 'quad') )

    n = size(dom.x{1}, 1);
    x = cat(3, dom.x{:});
    y = cat(3, dom.y{:});
    z = cat(3, dom.z{:});
    X = [ x(:,1,:) ; nan(1,1,ne) ; x(:,n,:) ; nan(1,1,ne) ; permute(x(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(x(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
    Y = [ y(:,1,:) ; nan(1,1,ne) ; y(:,n,:) ; nan(1,1,ne) ; permute(y(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(y(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
    Z = [ z(:,1,:) ; nan(1,1,ne) ; z(:,n,:) ; nan(1,1,ne) ; permute(z(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(z(n,:,:), [2 1 3]) ; nan(1,1,ne) ];

    if ( ~strcmpi(showEdges, 'off') )
        plot3(X(:), Y(:), Z(:), defaultStyle{:}, varargin{:})
    end

    % If the plot is not being added to another then plot the surface so
    % that the lines are more easily discernable.
    if ( (~holdState && strcmpi(showSurface, 'auto')) || strcmpi(showSurface, 'on') )
        % Plot the surface, making it slightly smaller so lines show up
        % more clearly.
        hold on
        for k = 1:ne
            surface(dom.x{k} - scl*vn{k}(:,:,1), ...
                    dom.y{k} - scl*vn{k}(:,:,2), ...
                    dom.z{k} - scl*vn{k}(:,:,3), ...
                    0*dom.x{k},                  ...
                    surfaceStyle{:});
        end
    end

end

if ( ~holdState )
    view(3)
    axis equal
    grid on
end

if ( ~holdState )
    hold off
end

end
