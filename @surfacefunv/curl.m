function g = curl(f)
%CURL   Curl of a SURFACEFUNV.
%   G = CURL(F) returns the surface curl of the SURFACEFUNV F as a
%   SURFACEFUNV. This operation only makes sense if F is tangent to the
%   surface.
%
% See also DIVERGENCE.

% Empty check:
if ( isempty(f) )
    g = surfacefunv;
    return
end

% Should we warn the user if f is not tangent to the surface?
fc = f.components;
g = surfacefunv( diff(fc{3}, 1, 2) - diff(fc{2}, 1, 3), ... % d/dy(fz) - d/dz(fy)
                 diff(fc{1}, 1, 3) - diff(fc{3}, 1, 1), ... % d/dz(fx) - d/dx(fz)
                 diff(fc{2}, 1, 1) - diff(fc{1}, 1, 2) );   % d/dx(fy) - d/dy(fx)

end
