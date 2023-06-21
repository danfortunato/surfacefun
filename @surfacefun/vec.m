function x = vec(f)
%VEC   Flatten a SURFACEFUN into a vector.
%   X = F.VEC() concatenates and flattens the function values from the
%   SURFACEFUN F to yield a column vector X of length NUMEL(F). If F is an
%   array-valued SURFACEFUN with M columns, then F.VEC() will be a matrix
%   of size NUMEL(F(1)) x M.
%
%   See also SUBSASGN.

n = order(f(1).domain) + 1;
numPatches = length(f(1).domain);
numDOF = numel(f(1));
numFuns = size(f, 2);

x = zeros(numDOF, numFuns);
for j = 1:numFuns
    for k = 1:numPatches
        x((k-1)*n^2+(1:n^2), j) = f(:,j).vals{k}(:);
    end
end

end
