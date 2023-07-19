function surf(f, varargin)
%SURF   3-D colored surface of a SURFACEFUN.
%   SURF(F) plots a colored parametric surface whose color is defined by
%   the values of the SURFACEFUN F.
%
%   SURF(..., 'PropertyName', PropertyValue, ...) sets the value of the
%   specified surface property. Multiple property values can be set with a
%   single statement.
%
%   See also PLOT, CONTOUR, MESH.

nf = builtin('numel', f);
if ( nf > 1 )
    % Plot each column of a surfacefun array in its own figure:
    for k = 1:nf
        figure(k)
        surf(f(k), varargin{:})
    end
    % Align the figures:
    alignfigs
    return
end

p = inputParser;
p.addParameter('plotpts', 100);
p.KeepUnmatched = true;
p.parse(varargin{:});
nplotpts = p.Results.plotpts;
argnames = fieldnames(p.Unmatched);
argvals = struct2cell(p.Unmatched);
args = [argnames(:).' ; argvals(:).'];
varargin = args(:).';

holdState = ishold();

for k = 1:length(f)
    u = f.vals{k};
    if ( ~isreal(u) )
        u = real(u);
    end
    x = f.domain.x{k};
    y = f.domain.y{k};
    z = f.domain.z{k};

    u = chebvals2plotvals(u, nplotpts);
    x = chebvals2plotvals(x, nplotpts);
    y = chebvals2plotvals(y, nplotpts);
    z = chebvals2plotvals(z, nplotpts);

    surface(x, y, z, u, varargin{:});
    hold on
end

if ( ~holdState )
    view(3)
    shading interp
    axis equal
    grid on
    hold off
end

end
