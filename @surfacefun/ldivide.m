function f = ldivide(f, g)
%.\   Pointwise left divide for SURFACEFUN.
%   F.\G divides G by F, where F and G may be SURFACEFUN objects or
%   scalars.
%
%   See also RDIVIDE, COMPOSE.

f = rdivide(g, f);

end
