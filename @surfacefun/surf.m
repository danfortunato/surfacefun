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

if ( ~all(f.domain.ptype == f.domain.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

switch ( f.domain.ptype(1) )

    case 'tri'

        x = f.domain.x;
        y = f.domain.y;
        z = f.domain.z;
        u = f.vals;

        for k = 1:length(f)
            x{k} = koornvals2plotvals(x{k}, nplotpts);
            y{k} = koornvals2plotvals(y{k}, nplotpts);
            z{k} = koornvals2plotvals(z{k}, nplotpts);
            u{k} = koornvals2plotvals(u{k}, nplotpts);
        end

        x_all = reshape([x{:}], [], 1);
        y_all = reshape([y{:}], [], 1);
        z_all = reshape([z{:}], [], 1);
        u_all = reshape([u{:}], [], 1);

        if ( ~isreal(u_all) )
            u_all = real(u_all);
        end

        T = trilattice(nplotpts);
        ntri = size(T, 1);
        ne = length(u);
        npts = length(u{1});
        T_all = repmat(T, [1 1 ne]) + reshape((0:ne-1)*npts, 1, 1, []);
        T_all = permute(T_all, [1 3 2]);
        T_all = reshape(T_all, ntri*ne, 3);

        hold on
        patch('Faces',    T_all,               ...
              'Vertices', [x_all y_all z_all], ...
              'FaceVertexCData', u_all,        ...
              varargin{:});

    case 'quad'

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
end

if ( ~holdState )
    view(3)
    shading interp
    %axis equal
    set(gca, 'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatioMode', 'auto')
    grid on
    hold off
end

end
