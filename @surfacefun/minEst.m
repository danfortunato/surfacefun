function m = minEst(f)
%MINEST   Estimate the minimum value of a SURFACEFUN.
%
%   See also MAXEST, NORM.

m = min(cellfun(@(u) min(u(:)), f.vals));

end
