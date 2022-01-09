function out = isempty(f)
%ISEMPTY   Test for empty SURFACEFUNV.
%   ISEMPTY(F) returns 1 if F is an empty SURFACEFUNV and 0 otherwise.

out = isempty(f.components) || ( isempty(f.components{1}) && ...
                                 isempty(f.components{2}) && ...
                                 isempty(f.components{3}) );

end
