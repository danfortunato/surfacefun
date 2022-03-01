function wireframe(dom, varargin)
%WIREFRAME   Plot the wireframe of a surface.

parser = inputParser;
parser.KeepUnmatched = true;
parser.addParameter('surface', 'auto', @(s) contains(lower(s), {'auto', 'on', 'off'}));
parse(parser, varargin{:});
showSurface = parser.Results.surface;
varargin = namedargs2cell(parser.Unmatched);

defaultStyle = {'Color', 'k', 'LineStyle', '-', 'LineWidth', 1};

holdState = ishold();

x = dom.x;
y = dom.y;
z = dom.z;
vn = dom.facenormals;

n = size(x{1}, 1);
ne = length(x);
x = cat(3, x{:});
y = cat(3, y{:});
z = cat(3, z{:});
X = [ x(:,1,:) ; nan(1,1,ne) ; x(:,n,:) ; nan(1,1,ne) ; permute(x(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(x(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
Y = [ y(:,1,:) ; nan(1,1,ne) ; y(:,n,:) ; nan(1,1,ne) ; permute(y(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(y(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
Z = [ z(:,1,:) ; nan(1,1,ne) ; z(:,n,:) ; nan(1,1,ne) ; permute(z(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(z(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
plot3(X(:), Y(:), Z(:), defaultStyle{:}, varargin{:})

% If the plot is not being added to another then plot the surface so
% that the lines are more easily discernable.
if ( (~holdState && strcmpi(showSurface, 'auto')) || strcmpi(showSurface, 'on') )
    % Plot the surface, making it slightly smaller so lines show up
    % more clearly.
    hold on
    for k = 1:length(dom)
        scl = 0.01;
        surface(dom.x{k} - scl*vn{k}(:,:,1), ...
                dom.y{k} - scl*vn{k}(:,:,2), ...
                dom.z{k} - scl*vn{k}(:,:,3), ...
                0*dom.x{k}, 'FaceColor', 'w', 'EdgeColor', 'k', ...
                'AmbientStrength', 0.6, 'DiffuseStrength', 0.4, 'SpecularStrength', 0.3);
    end
end

if ( ~holdState )
    axis equal
end

if ( ~holdState )
    hold off
end

end
