function mesh(dom, varargin)
%MESH   Plot the mesh of a surface.

holdState = ishold();

x = dom.x;
y = dom.y;
z = dom.z;
vn = dom.facenormals;

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
                0*x{k}, 'FaceColor', 'w', 'EdgeColor', 'None');
    end
end

if ( ~holdState )
    axis equal
    hold off
end

end
