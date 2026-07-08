function f = compose(op, f, g)
%COMPOSE   Compose a function handle with a SURFACEFUN.
%   COMPOSE(OP, F) returns a SURFACEFUN representing OP(F), where OP is a
%   function handle and F is a SURFACEFUN. The function handle OP is
%   applied to F in value space.
%
%   COMPOSE(OP, F, G) returns a SURFACEFUN representing OP(F, G), where OP
%   is a function handle and at least one of F and G is a SURFACEFUN. The
%   function handle OP is applied to F and G in value space.

if ( nargin == 2 )
    nf = builtin('numel', f);
    for k = 1:nf
       f(k).vals = cellfun(op, f(k).vals, 'UniformOutput', false);
    end
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefun') )
    % F and G are both SURFACEFUNs:
    nf = builtin('numel', f);
    if ( ~all(size(f) == size(g)) )
        error('SURFACEFUN:compose:dims', 'Matrix dimensions must agree.');
    end
    for k = 1:nf
        f(k).vals = cellfun(op, f(k).vals, g(k).vals, 'UniformOutput', false);
    end
elseif ( isa(f, 'surfacefun') )
    % G is not a SURFACEFUN:
    nf = builtin('numel', f);
    if ( isscalar(g) )
        for k = 1:nf
            f(k).vals = cellfun(@(x) op(x,g), f(k).vals, 'UniformOutput', false);
        end
    elseif ( all(size(g) == size(f)) )
        for k = 1:nf
            f(k).vals = cellfun(@(x) op(x,g(k)), f(k).vals, 'UniformOutput', false);
        end
    else
        error('SURFACEFUN:compose:dims', 'Matrix dimensions must agree.');
    end
elseif ( isa(g, 'surfacefun') )
    % F is not a SURFACEFUN:
    ng = builtin('numel', g);
    if ( isscalar(f) )
        for k = 1:ng
            g(k).vals = cellfun(@(x) op(f,x), g(k).vals, 'UniformOutput', false);
        end
    elseif ( all(size(f) == size(g)) )
        for k = 1:ng
            g(k).vals = cellfun(@(x) op(f(k),x), g(k).vals, 'UniformOutput', false);
        end
    else
        error('SURFACEFUN:compose:dims', 'Matrix dimensions must agree.');
    end
    f = g;
end

end
