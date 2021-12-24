function contour(f, varargin)
%CONTOUR   Contour plot of a SURFACEFUN.
%   CONTOUR(F) is a contour plot of F treating the values of F as heights
%   above or below the surface on which F is defined. Contours are the
%   level curves of F for some values V. The values V are chosen
%   automatically.
%
%   CONTOUR(F, N) draws N contour lines, choosing the levels automatically.
%   
%   CONTOUR(F, V) draws a contour line for each level specified in the
%   vector V. Use CONTOUR(F, [V V]) to compute a single contour at the
%   level V.
%
%   See also PLOT, SURF.

holdState = ishold();
N = 10;
levels = [];

% See if an N was given:
if ( nargin > 1 && isnumeric(varargin{1}) )
    v1 = varargin{1};
    if ( isscalar(v1) )
        N = v1;
    else
        levels = v1;
    end
    varargin(1) = [];
end

% Determine some levels:
if ( isempty(levels) )
    minu = minEst(f);
    maxu = maxEst(f);
    levels = linspace(minu, maxu, N);
end

% Loop over the patches:
m = 100;
[uu, vv] = meshgrid(linspace(-1, 1, m));
for j = 1:length(f)
    u = chebvals2plotvals(f.vals{j});
    if ( ~isreal(u) )
        u = abs(u);
    end
    
    % Get contour lines.
    [C, H] = contour(uu, vv, u, levels, varargin{:});
    
    % Extract out the options we need to plot the contours with plot3.
    lw = H.LineWidth;
    ls = H.LineStyle;
    lc = H.LineColor;
    levelList = H.LevelList;
    clrmap = parula(numel(levelList));

    % Remove the contour plot that was generated.
    delete(H);

    % If the plot is not being added to another then plot the surface so
    % that the lines are more easily discernable.
    if ( ~holdState )
        % Plot the surface, making it slightly smaller so lines show up
        % more clearly.
        xx = f.domain.x{j};
        yy = f.domain.y{j};
        zz = f.domain.z{j};
        scl = 0.99;
        surf(scl*xx, scl*yy, scl*zz, 1+0*xx, 'FaceColor', 'w', 'EdgeColor', 'None');
        hold on
    end

    % Plot the contours on the surface.
    k = 1;
    while ( k < size(C, 2) )
        kl = C(2, k);
        v = k+1:k+kl;
        xv = bary2d(f.domain.x{j}, C(1, v), C(2, v));
        yv = bary2d(f.domain.y{j}, C(1, v), C(2, v));
        zv = bary2d(f.domain.z{j}, C(1, v), C(2, v));

        % If the line color is a float then we are plotting all contours in
        % a single color.
        if ( isfloat(lc) )
            plot3(xv, yv, zv, 'LineWidth', lw, 'Color', lc, 'LineStyle', ls);
        else
            % We need to plot each contour in a color using the default 
            % colormap. Determine the color for the level being plotted.
            clr = clrmap(abs(C(1, k) - levelList) < 10*eps, :);
            plot3(xv, yv, zv, 'LineWidth', lw, 'Color', clr, 'LineStyle', ls);
        end        
        k = k+kl+1;
        hold on
    end
end
axis equal

if ( ~holdState )
    hold off
end

end

function out = bary2d(vals, x, y)

out = zeros(size(x));
yvals = bary(y(:), vals).';
for k = 1:numel(x)
    out(k) = bary(x(k), yvals(:,k));
end
out = reshape(out, size(x));

end
