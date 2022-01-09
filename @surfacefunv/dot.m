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

fc = f.components;
gc = g.components;
h = fc{1}.*gc{1} + fc{2}.*gc{2} + fc{3}.*gc{3};

end
