function f = rdivide(f, g)
%./   Pointwise right divide for SURFACEFUN.
%   F./G divides F by G, where F and G may be SURFACEFUN objects or
%   scalars.
%
%   See also LDIVIDE, COMPOSE.

% If either f or g are empty then return an empty SURFACEFUN object.
if ( isempty(f) )
    return
elseif ( isempty(g) )
    f = g;
    return
end

if ( isa(f, 'surfacefun') && isa(g, 'surfacefun') )
    % Divide two SURFACEFUNs:
    % TODO: Check that F and G have the same domain.
    [ss, wzero] = singleSignTest(g);
    if ( ss && ~wzero )
        f = compose(@rdivide, f, g);
    else
        error('SURFACEFUN:rdivide:zero', ...
              'Attempting to invert a surfacefun with a root.');
    end
elseif ( isa(f, 'surfacefun') && isnumeric(g) && isscalar(g) )
    % Divide SURFACEFUN F by scalar G:
    f.vals = cellfun(@(vals) (1./g).*vals, f.vals, 'UniformOutput', false);
elseif ( isnumeric(f) && isscalar(f) && isa(g, 'surfacefun') )
    % Divide scalar F by SURFACEFUN G:
    [ss, wzero] = singleSignTest(g);
    if ( ss && ~wzero )
        f = compose(@rdivide, f, g);
    else
        error('SURFACEFUN:rdivide:zero', ...
              'Attempting to invert a surfacefun with a root.');
    end
else
    error('SURFACEFUN:rdivide:invalid', ...
        'F and G must be scalars or surfacefun objects.')
end

end
