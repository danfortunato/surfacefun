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
    ff = cellfun(op, f.vals, 'UniformOutput', false);
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefun') )
    % F and G are both SURFACEFUNs:
    ff = cellfun(op, f.vals, g.vals, 'UniformOutput', false);
elseif ( isa(f, 'surfacefun') )
    % G is not a SURFACEFUN:
    ff = cellfun(@(x) op(x,g), f.vals, 'UniformOutput', false);
elseif ( isa(g, 'surfacefun') )
    % F is not a SURFACEFUN:
    ff = cellfun(@(x) op(f,x), g.vals, 'UniformOutput', false);
    f = g;
end

f.vals = ff;

end
