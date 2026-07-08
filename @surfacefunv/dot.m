function h = dot(f, g)
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of the SURFACEFUNV objects F and G.
%   DOT(F, G) is the same as F'*G.
%
% See also CROSS.

if ( isempty(f) || isempty(g) )
    h = surfacefun;
    return
end

if ( isnumeric(f) )
    h = dot(g, f);
    return
elseif ( isnumeric(g) )
    if ( numel(g) == 3 )
        fc = f.components;
        h = fc{1}.*g(1) + fc{2}.*g(2) + fc{3}.*g(3);
    else
        error('SURFACEFUNV:dot:invalid', ...
            'F and G must be surfacefunv objects or constant vectors.');
    end
elseif ( isa(f, 'surfacefunv') && isa(g, 'surfacefunv') )
    fc = f.components;
    gc = g.components;
    h = fc{1}.*gc{1} + fc{2}.*gc{2} + fc{3}.*gc{3};
else
    error('SURFACEFUNV:dot:invalid', ...
        'F and G must be surfacefunv objects or constant vectors');
end

end
