function out = isempty(f)
%ISEMPTY   Test for empty SURFACEFUN.
%   ISEMPTY(F) returns 1 if F is an empty SURFACEFUN and 0 otherwise.

out = isempty(f.vals);

end
