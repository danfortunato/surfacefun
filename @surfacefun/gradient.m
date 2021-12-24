function [fx, fy, fz] = gradient(f)
%GRADIENT   Gradient of a SURFACEFUN.
%   [FX, FY, FZ] = GRADIENT(F) returns the numerical gradient of the
%   SURFACEFUN F, where FX, FY, and FZ are the derivatives of F in the x-,
%   y-, and z-directions. The derivatives are returned as SURFACEFUN
%   objects.
%
%   See also GRAD.

fx = diff(f, 1, 1);
fy = diff(f, 1, 2);
fz = diff(f, 1, 3);

end
