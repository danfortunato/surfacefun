function f = times(f, g)
%.*   Pointwise multiplication for SURFACEFUN.
%   F.*G multiplies F and G pointwise, where F and G may be SURFACEFUN
%   objects or scalars.
%
%   See also MTIMES, COMPOSE.

if ( isa(g, 'surfacefunv') )
    h = g;
    h.components{1} = times(f, g.components{1});
    h.components{2} = times(f, g.components{2});
    h.components{3} = times(f, g.components{3});
    f = h;
    return
end

if ( ~isa(f, 'surfacefun') )
    % Ensure F is the SURFACEFUN:
    f = times(g, f);
    return
elseif ( isa(g, 'surfacefun' ) )
    % Multiply two SURFACEFUNs:
    % TODO: Check that F and G have the same domain.
    if ( ~all(size(f) == size(g)) )
        error('SURFACEFUN:times:dims', 'Matrix dimension must agree.');
    end
    nfuns = prod(size(f)); %#ok<PSIZE>
    for k = 1:nfuns
        f(k) = compose(@times, f(k), g(k));
    end
elseif ( isnumeric(g) && isscalar(g) )
    % Multiply SURFACEFUN F by scalar G:
    f.vals = cellfun(@(vals) g*vals, f.vals, 'UniformOutput', false);
elseif ( isnumeric(g) )
    % Multiply SURFACEFUN F by matrix G:
    if ( ~all(size(f) == size(g)) )
        error('SURFACEFUN:times:dims', 'Matrix dimension must agree.');
    end
    nfuns = prod(size(f)); %#ok<PSIZE>
    for k = 1:nfuns
        f(k) = g(k)*f(k);
    end
else
    error('SURFACEFUN:times:invalid', ...
        'F and G must be scalars or surfacefun objects.')
end

end
