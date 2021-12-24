function N = numel(f)
%NUMEL   Number of degrees of freedom in a SURFACEFUN.
%   N = NUMEL(F) returns the total number of degrees of freedom in the
%   SURFACEFUN object F.

N = 0;
for k = 1:length(f)
    N = N + numel(f.vals{k});
end

end
