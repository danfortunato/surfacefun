function f = sqrt(f)
%SQRT   Square root of a SURFACEFUN.
%   SQRT(F) returns the square root of the SURFACEFUN F.
%
%   See also POWER, COMPOSE.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@sqrt, f);

end
