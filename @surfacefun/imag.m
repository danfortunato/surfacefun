function sol = imag(sol)
%IMAG   Imaginary part of a SURFACEFUN.
%
%   See also REAL.

sol = compose(@imag, sol);

end
