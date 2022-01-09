function varargout = gradient(f)
%GRADIENT   Gradient of a SURFACEFUN.
%   [FX, FY, FZ] = GRADIENT(F) returns the numerical gradient of the
%   SURFACEFUN F, where FX, FY, and FZ are the derivatives of F in the x-,
%   y-, and z-directions. The derivatives are returned as SURFACEFUN
%   objects unless one output is requested, in which case the gradient is
%   returned as a SURFACEFUNV.
%
%   See also GRAD.

fx = diff(f, 1, 1);
fy = diff(f, 1, 2);
fz = diff(f, 1, 3);

if ( nargout == 0 || nargout == 1 )
    varargout{1} = surfacefunv(fx, fy, fz);
else
    varargout = {fx, fy, fz};
end

end
