function h = times(f, g)
%.*   Pointwise multiply for SURFACEFUN(V).
%   F.*G divides F by G, where F and G may be SURFACEFUN(V) objects or
%   scalars.

% Empty check:
h = surfacefunv;
if ( isempty(f) || isempty(g) )
    return
end

if ( isa(f, 'surfacefunv') && isa(g, 'surfacefunv') )
    % Multiply two SURFACEFUNVs:
    h.components{1} = f.components{1} .* g.components{1};
    h.components{2} = f.components{2} .* g.components{2};
    h.components{3} = f.components{3} .* g.components{3};
elseif ( isa(f, 'surfacefunv') && isa(g, 'surfacefun') )
    % Multiply SURFACEFUNV F by SURFACEFUN G:
    h.components{1} = f.components{1} .* g;
    h.components{2} = f.components{2} .* g;
    h.components{3} = f.components{3} .* g;
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefunv') )
    % Multiply SURFACEFUN F by SURFACEFUNV G:
    h.components{1} = f .* g.components{1};
    h.components{2} = f .* g.components{2};
    h.components{3} = f .* g.components{3};
elseif ( isa(f, 'surfacefunv') && isnumeric(g) )
    if ( isscalar(g) )
        % Multiply SURFACEFUNV F by scalar G:
        h.components{1} = f.components{1} .* g;
        h.components{2} = f.components{2} .* g;
        h.components{3} = f.components{3} .* g;
    elseif ( numel(g) == 3 )
        % Multiply SURFACEFUNV F by vector G:
        h.components{1} = f.components{1} .* g(1);
        h.components{2} = f.components{2} .* g(2);
        h.components{3} = f.components{3} .* g(3);
    else
        error('SURFACEFUNV:times:invalid', ...
            'F and G must be surfacefunv objects, scalars, or constant vectors.');
    end
elseif ( isnumeric(f) && isa(g, 'surfacefunv') )
    if ( isscalar(f) )
        % Multiply scalar F by SURFACEFUNV G:
        h.components{1} = f .* g.components{1};
        h.components{2} = f .* g.components{2};
        h.components{3} = f .* g.components{3};
    elseif ( numel(f) == 3 )
        % Multiply vector F by SURFACEFUNV G:
        h.components{1} = f(1) .* g.components{1};
        h.components{2} = f(2) .* g.components{2};
        h.components{3} = f(3) .* g.components{3};
    else
        error('SURFACEFUNV:times:invalid', ...
            'F and G must be surfacefunv objects, scalars, or constant vectors.');
    end
else
    error('SURFACEFUNV:times:invalid', ...
        'F and G must be surfacefunv objects, scalars, or constant vectors.');
end

end

