function m = mean2(f)
%MEAN2   Mean of a SURFACEFUN.
%   MEAN2(F) returns the mean of the SURFACEFUN F.

m = integral2(f) / surfacearea(f.domain);

end
