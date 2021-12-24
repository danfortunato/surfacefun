function m = maxEst(f)
%MAXEST   Estimate the maximum value of a SURFACEFUN.
%
%   See also MINEST, NORM.

m = max(cellfun(@(u) max(u(:)), f.vals));

end
