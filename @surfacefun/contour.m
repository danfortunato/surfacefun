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

ff = f.vals;
xx = f.domain.x;
yy = f.domain.y;
zz = f.domain.z;
n = order(f)+1;
[bk, ~, vk] = chebpts(n);

if ( ~all(f.domain.ptype == f.domain.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

if ( f.domain.ptype(1) == 'tri' ) %#ok<BDSCA>
    % Convert triangle points to tensor-product points using the Duffy
    % transformation so that we can use MATLAB's built-in contour().
    [eta1, eta2] = chebpts2(n, n, [0 1 0 1]);
    xd = eta1.*(1-eta2);
    yd = eta2;
    K  = koornwinder(n-1);
    Kd = koornwinder(n-1, xd, yd);
    B = Kd / K;
    for k = 1:length(f)
        ff{k} = reshape(B*ff{k}, n, n);
        xx{k} = reshape(B*xx{k}, n, n);
        yy{k} = reshape(B*yy{k}, n, n);
        zz{k} = reshape(B*zz{k}, n, n);
    end
end

% Loop over the patches:
m = 100;
[uu, vv] = meshgrid(linspace(-1, 1, m));
for j = 1:length(f)
    u = chebvals2plotvals(ff{j});
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
    %if ( ~holdState )
        % Plot the surface, making it slightly smaller so lines show up
        % more clearly.
        %xx = f.domain.x{j};
        %yy = f.domain.y{j};
        %zz = f.domain.z{j};
        %scl = 0.99;
        %surf(scl*xx{j}, scl*yy{j}, scl*zz{j}, 1+0*xx{j}, 'FaceColor', 'w', 'EdgeColor', 'None');
        %hold on
    %end

    % Plot the contours on the surface.
    k = 1;
    while ( k < size(C, 2) )
        kl = C(2, k);
        v = k+1:k+kl;
        xv = bary2d(xx{j}, C(1, v), C(2, v), bk, vk);
        yv = bary2d(yy{j}, C(1, v), C(2, v), bk, vk);
        zv = bary2d(zz{j}, C(1, v), C(2, v), bk, vk);

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

% Plot the surface, making it slightly smaller so lines show up
plot(f.domain, edges='off', surface='on')

axis equal

if ( ~holdState )
    hold off
end

end

function out = bary2d(vals, x, y, bk, vk)

out = zeros(size(x));
yvals = bary(y(:), vals, bk, vk).';
for k = 1:numel(x)
    out(k) = bary(x(k), yvals(:,k), bk, vk);
end
out = reshape(out, size(x));

end
