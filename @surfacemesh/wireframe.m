function wireframe(dom, varargin)
%WIREFRAME   Plot the wireframe of a surface.

holdState = ishold();

x = dom.x;
y = dom.y;
z = dom.z;

n = size(x{1}, 1);
ne = length(x);
x = cat(3, x{:});
y = cat(3, y{:});
z = cat(3, z{:});
X = [ x(:,1,:) ; nan(1,1,ne) ; x(:,n,:) ; nan(1,1,ne) ; permute(x(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(x(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
Y = [ y(:,1,:) ; nan(1,1,ne) ; y(:,n,:) ; nan(1,1,ne) ; permute(y(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(y(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
Z = [ z(:,1,:) ; nan(1,1,ne) ; z(:,n,:) ; nan(1,1,ne) ; permute(z(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(z(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
plot3(X(:), Y(:), Z(:), 'k-', 'LineWidth', 1, varargin{:})

% If the plot is not being added to another then plot the surface so
% that the lines are more easily discernable.
if ( ~holdState )
    % Plot the surface, making it slightly smaller so lines show up
    % more clearly.
    hold on
    for k = 1:length(dom)
        scl = 0.99;
        surf(scl*dom.x{k}, scl*dom.y{k}, scl*dom.z{k}, 1+0*dom.x{k}, 'FaceColor', 'w', 'EdgeColor', 'None');
    end
    axis equal
end

if ( ~holdState )
    hold off
end

end
