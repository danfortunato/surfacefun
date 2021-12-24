function h = plus(f, g)
%+   Plus for SURFACEFUN.
%   F + G adds the SURFACEFUN F and G. F and G must have the same domains
%   and discretization sizes. F and G may also be scalars.
%
%   See also MINUS.

if ( isnumeric( f ) )
    h = plus(g, f);
    return
elseif ( isnumeric( g ) )
    h = f;
    for k = 1:size(h.vals,1)
        h.vals{k} = h.vals{k} + g;
    end
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefun') )
    h = f;
    % TODO: Assume on the same grid for now.
    h.vals = cellfun(@plus, f.vals , g.vals, 'UniformOutput', false);
end

end
