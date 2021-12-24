function A = surfacearea(dom)
%SURFACEAREA   Compute the surface area of a surface.
%
%   See also VOLUME.

if ( isempty(dom) )
    A = 0;
    return
end

I = surfacefun(@(x,y,z) 1+0*x, dom);
A = integral2(I);

end
