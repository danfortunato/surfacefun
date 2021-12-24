function f = mtimes(f, c)
%*   Scale a SURFACEFUN.
%   c*F or F*c multiplies a SURFACEFUN F by a scalar c.
%
%   See also TIMES.

if ( ~isa(f, 'surfacefun') )
    % Ensure F is the SURFACEFUN:
    f = mtimes(c, f);
    return
elseif ( isa(c, 'surfacefun' ) )
    % MTIMES should not be used to multiply two SURFACEFUNs:
    error('SURFACEFUN:mtimes:twofuns', ...
        ['Cannot multiply two surfacefuns with ''*''. ', ...
         'Did you mean ''.*''?\n'])
elseif ( isnumeric(c) && isscalar(c) )
    % Multiply SURFACEFUN F by scalar c:
    f.vals = cellfun(@(vals) c*vals, f.vals, 'UniformOutput', false);
else
    error('SURFACEFUN:mtimes:invalid', 'c must be a scalar.')
end

end
