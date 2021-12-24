function I = integral(f, varargin)
%INTEGRAL   Double integral of a SURFACEFUN.
%   I = INTEGRAL(F) returns the double integral of the SURFACEFUN F over
%   its domain.
%
%   I = INTEGRAL(F, 'all') returns an array of double integrals over each
%   patch of F.

I = integral2(f, varargin);

end
