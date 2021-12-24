function I = sum2(f, varargin)
%SUM2   Double integral of a SURFACEFUN.
%   I = SUM2(F) returns the double integral of the SURFACEFUN F over the
%   surface.
%
%   I = SUM2(F, 'all') returns an array of double integrals over each patch
%   of F.
%
%   See also INTEGRAL2.

I = integral2(f, varargin);

end
