function h = plus(f, g)
%+   Plus for SURFACEFUNV.
%   F + G adds the SURFACEFUNV F and G. F and G must have the same domains
%   and discretization sizes. F and G may also be scalars.
%
%   See also MINUS.

if ( isnumeric(f) )
    h = plus(g, f);
    return
elseif ( isnumeric(g) )
    h = f;
    if ( numel(g) == 1 )
        h.components{1} = f.components{1} + g;
        h.components{2} = f.components{2} + g;
        h.components{3} = f.components{3} + g;
    elseif ( numel(g) == 3 )
        h.components{1} = f.components{1} + g(1);
        h.components{2} = f.components{2} + g(2);
        h.components{3} = f.components{3} + g(3);
    else
        error('SURFACEFUNV:plus:invalid', ...
            'F and G must be surfacefunv objects, constant vectors, or scalars.');
    end
elseif ( isa(f, 'surfacefunv') && isa(g, 'surfacefunv') )
    h = f;
    h.components{1} = f.components{1} + g.components{1};
    h.components{2} = f.components{2} + g.components{2};
    h.components{3} = f.components{3} + g.components{3};
else
    error('SURFACEFUNV:plus:invalid', ...
        'F and G must be surfacefunv objects, constant vectors, or scalars.');
end

end
