function d = div(f)
%DIV   Surface divergence of a SURFACEFUNV.
%   D = DIV(F) returns the surface divergence of the SURFACEFUNV F as a
%   SURFACEFUN. This operation only makes sense if F is tangent to the
%   surface.
%
%   This is shorthand for DIVERGENCE(F).
%
% See also DIVERGENCE, CURL.

d = divergence(f);

end
