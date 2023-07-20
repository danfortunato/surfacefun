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
    f = compose(@times, f, g);
elseif ( isnumeric(g) && isscalar(g) )
    % Multiply SURFACEFUN F by scalar G:
    f = compose(@(x) g*x, f);
elseif ( isnumeric(g) )
    % Multiply SURFACEFUN F by matrix G:
    if ( ~all(size(f) == size(g)) )
        error('SURFACEFUN:times:dims', 'Matrix dimensions must agree.');
    end
    nf = builtin('numel', f);
    for k = 1:nf
        f(k) = g(k)*f(k);
    end
else
    error('SURFACEFUN:times:invalid', ...
        'F and G must be scalars or surfacefun objects.')
end

end
