function f = mrdivide(f, c)
%/   Right scalar divide for SURFACEFUN.
%   F/C divides the SURFACEFUN F by a scalar C.
%
%   See also RDIVIDE, MLDIVIDE.

if ( isnumeric(c) )
    f = rdivide(f, c);
else
    error('SURFACEFUN:mrdivide:mrdivide', 'Not supported. Did you mean ./ ?');
end

end
