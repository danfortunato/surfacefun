function [fx, fy, fz] = grad(f)
%GRAD   Gradient of a SURFACEFUN.
%   [FX, FY] = GRAD(F) returns the numerical gradient of the SURFACEFUN F,
%   where FX, FY, and FZ are the derivatives of F in the x-, y-, and
%   z-directions. The derivatives are returned as SURFACEFUN objects.
%
%   This is shorthand for GRADIENT(F).
%
%   See also GRADIENT.

[fx, fy, fz] = gradient(f);

end
