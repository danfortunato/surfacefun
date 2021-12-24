function sol = real(sol)
%REAL   Real part of a SURFACEFUN.
%
%   See also IMAG.

sol = compose(@real, sol);

end
