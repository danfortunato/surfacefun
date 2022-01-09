function varargout = grad(f)
%GRAD   Gradient of a SURFACEFUN.
%   [FX, FY, FZ] = GRAD(F) returns the numerical gradient of the SURFACEFUN
%   F, where FX, FY, and FZ are the derivatives of F in the x-, y-, and
%   z-directions. The derivatives are returned as SURFACEFUN objects
%   unless one output is requested, in which case the gradient is returned
%   as a SURFACEFUNV.
%
%   This is shorthand for GRADIENT(F).
%
%   See also GRADIENT.

[varargout{1:nargout}] = gradient(f);

end
