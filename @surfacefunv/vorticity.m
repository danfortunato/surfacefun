function v = vorticity(f)
%VORTICITY   Surface vorticity of a SURFACEFUNV.
%   V = VORTICITY(F) returns the surface vorticity of the SURFACEFUNV F as
%   a SURFACEFUN. Vorticity is the defined as the normal component of
%   CURL(F) on the surface. It is the generalization of the standard 2D
%   scalar vorticity to a surface.
%
% See also VORT, DIV, CURL.

% Empty check:
if ( isempty(f) )
    v = surfacefun;
    return
end

% Vorticity is n dot curl(F), where n is the unit normal to the surface.
v = dot(normal(domain(f)), curl(f));

end
