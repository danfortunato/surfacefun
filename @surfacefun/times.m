function f = times(f, g)
%.*   Pointwise multiplication for SURFACEFUN.
%   F.*G multiplies F and G, where F and G may be SURFACEFUN objects or
%   scalars.
%
%   See also MTIMES, COMPOSE.

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
    f.vals = cellfun(@(vals) g*vals, f.vals, 'UniformOutput', false);
else
    error('SURFACEFUN:times:invalid', ...
        'F and G must be scalars or surfacefun objects.')
end

end
