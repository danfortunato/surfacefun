function f = exp(f)
%EXP   Exponential of a SURFACEFUN.
%   EXP(F) returns the exponential of the SURFACEFUN F.
%
%   See also LOG.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@exp, f);

end
