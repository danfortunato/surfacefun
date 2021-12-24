function varargout = surf(f, varargin)
%SURF   3-D colored surface of a SURFACEFUN.
%   SURF(F) plots a colored parametric surface whose color is defined by
%   the values of the SURFACEFUN F.
%
%   SURF(..., 'PropertyName', PropertyValue, ...) sets the value of the
%   specified surface property. Multiple property values can be set with a
%   single statement.
%
%   H = SURF(...) returns a handle to a surface plot object.
%
%   See also PLOT, CONTOUR, MESH.

holdState = ishold();

for k = 1:length(f)
    u = f.vals{k};
    if ( ~isreal(u) )
        u = abs(u);
    end
    x = f.domain.x{k};
    y = f.domain.y{k};
    z = f.domain.z{k};
    h(k) = surf(x, y, z, u, 'EdgeAlpha', 1, varargin{:}); %#ok<AGROW>
    hold on
end

if ( ~holdState )
    shading interp
    axis equal
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
