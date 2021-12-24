function f = log(f)
%LOG   Natural logarithm of a SURFACEFUN.
%   LOG(F) returns the natural logarithm of the SURFACEFUN F. If F has
%   any roots in its domain, then the representation is likely to be
%   inaccurate.
%
%   See also LOG2, LOG10, EXP.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@log, f);

end
