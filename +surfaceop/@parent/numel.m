function N = numel(P)
%NUMEL   Number of degrees of freedom in a PARENT.
%   N = NUMEL(P) returns the total number of degrees of freedom in the
%   PARENT object P.

N = numel(P.child1) + numel(P.child2);

end
