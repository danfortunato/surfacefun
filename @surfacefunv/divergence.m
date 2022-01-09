function d = divergence(f)
%DIVERGENCE   Surface divergence of a SURFACEFUNV.
%   D = DIVERGENCE(F) returns the surface divergence of the SURFACEFUNV F
%   as a SURFACEFUN. This operation only makes sense if F is tangent to the
%   surface.
%
% See also CURL.

% Empty check:
if ( isempty(f) )
    d = surfacefun;
    return
end

% Should we warn the user if f is not tangent to the surface?
fc = f.components;
d = diff(fc{1}, 1, 1) + diff(fc{2}, 1, 2) + diff(fc{3}, 1, 3);

end
