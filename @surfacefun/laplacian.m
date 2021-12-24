function L = laplacian(f)
%LAPLACIAN   Laplacian of a SURFACEFUN.
%   LAPLACIAN(F) returns a SURFACEFUN representing the Laplacian of F.
%
%   See also LAP.

L = diff(f, 2, 2) + diff(f, 2, 1);

end
