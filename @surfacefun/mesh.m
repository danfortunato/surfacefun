function varargout = mesh(f, varargin)
%MESH   3-D mesh surface of a SURFACEFUN.
%   MESH(F) plots the tensor-product Chebyshev grids for the patches of F,
%   colored according to the values of F.
%
%   See also PLOT, SURF, CONTOUR.

holdState = ishold();

[x, y, z] = getGrid(f);
for k = 1:length(f)
    u = f.vals{k};
    if ( ~isreal(u) )
        u = abs(u);
    end
    h(k) = mesh(x{k}, y{k}, z{k}, u, varargin{:}); %#ok<AGROW>
    hold on
end

if ( ~holdState )
    axis equal
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
