function h = rdivide(f, g)
%./   Pointwise right divide for SURFACEFUN.
%   F./G divides F by G, where F and G may be SURFACEFUN objects or
%   scalars.
%
%   See also LDIVIDE, COMPOSE.

% Empty check:
h = surfacefunv;
if ( isempty(f) || isempty(g) )
    return
end

if ( isa(f, 'surfacefunv') && isa(g, 'surfacefunv') )
    % Divide two SURFACEFUNVs:
    h.components{1} = f.components{1} ./ g.components{1};
    h.components{2} = f.components{2} ./ g.components{2};
    h.components{3} = f.components{3} ./ g.components{3};
elseif ( isa(f, 'surfacefunv') && isa(g, 'surfacefun') )
    % Divide SURFACEFUNV F by SURFACEFUN G:
    h.components{1} = f.components{1} ./ g;
    h.components{2} = f.components{2} ./ g;
    h.components{3} = f.components{3} ./ g;
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefunv') )
    % Divide SURFACEFUN F by SURFACEFUNV G:
    h.components{1} = f ./ g.components{1};
    h.components{2} = f ./ g.components{2};
    h.components{3} = f ./ g.components{3};
elseif ( isa(f, 'surfacefunv') && isnumeric(g) )
    if ( isscalar(g) )
        % Divide SURFACEFUNV F by scalar G:
        h.components{1} = f.components{1} ./ g;
        h.components{2} = f.components{2} ./ g;
        h.components{3} = f.components{3} ./ g;
    elseif ( numel(g) == 3 )
        % Divide SURFACEFUNV F by vector G:
        h.components{1} = f.components{1} ./ g(1);
        h.components{2} = f.components{2} ./ g(2);
        h.components{3} = f.components{3} ./ g(3);
    else
        error('SURFACEFUNV:rdivide:invalid', ...
            'F and G must be surfacefunv objects, scalars, or constant vectors.');
    end
elseif ( isnumeric(f) && isa(g, 'surfacefunv') )
    if ( isscalar(f) )
        % Divide scalar F by SURFACEFUNV G:
        g.components{1} = f ./ g.components{1};
        g.components{2} = f ./ g.components{2};
        g.components{3} = f ./ g.components{3};
    elseif ( numel(f) == 3 )
        % Divide vector F by SURFACEFUNV G:
        g.components{1} = f(1) ./ g.components{1};
        g.components{2} = f(2) ./ g.components{2};
        g.components{3} = f(3) ./ g.components{3};
    else
        error('SURFACEFUNV:rdivide:invalid', ...
            'F and G must be surfacefunv objects, scalars, or constant vectors.');
    end
else
    error('SURFACEFUNV:rdivide:invalid', ...
        'F and G must be surfacefunv objects, scalars, or constant vectors.');
end

end
