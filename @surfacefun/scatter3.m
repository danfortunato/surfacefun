function varargout = scatter3(f, varargin)
%SCATTER3   3-D scatter plot of a SURFACEFUN.
%   SCATTER3(F) displays markers at the tensor-product Chebyshev grids for
%   the patches of F, colored according to the values of F.
%
%   H = SCATTER3(...) returns handles to the scatter objects created.
%
%   See also PLOT, MESH.

holdState = ishold();

[x, y, z] = getGrid(f);
for k = 1:length(f)
    h(k) = scatter3(x{k}(:), y{k}(:), z{k}(:), [], f.vals{k}(:), varargin{:}); %#ok<AGROW>
    hold on
end
axis equal

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
