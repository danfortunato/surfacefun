function f = power(f, g)
%.^   Pointwise power of a SURFACEFUN.
%   F.^G returns a SURFACEFUN F to the scalar power G, a scalar F to the
%   SURFACEFUN power G, or a SURFACEFUN F to the SURFACEFUN power G.
%
%   See also COMPOSE.

if ( isempty(f) )
    return
elseif ( isempty(g) )
    f = g;
    return
else
    f = compose(@power, f, g);
end

end
