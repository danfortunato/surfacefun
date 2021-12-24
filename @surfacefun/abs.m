function f = abs(f)
%ABS   Absolute value of a SURFACEFUN
%   ABS(F) returns the absolute value of the SURFACEFUN F.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@abs, f);

end
