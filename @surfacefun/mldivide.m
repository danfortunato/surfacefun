function f = mldivide(c, f)
%\   Left scalar divide for SURFACEFUN.
%   C\F divides the SURFACEFUN F by a scalar C.
%
%   See also LDIVIDE, MRDIVIDE.

if ( isnumeric(c) )
    f = ldivide(c, f);
else
    error('SURFACEFUN:mldivide:mldivide', 'Not supported. Did you mean .\ ?');
end

end
