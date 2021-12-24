function [x, y, z] = getGrid(f, kk)
%GETGRID   Get the grid of a SURFACEFUN.
%   [X, Y, Z] = GETGRID(F) returns the tensor-product Chebyshev grids (X,
%   Y, Z) for the patches of F.
%
%   [X, Y, Z] = GETGRID(F, KK) returns the tensor-product Chebyshev grids
%   for the KK-th patches of F.

if ( nargin == 1 )
    kk = 1:length(f);
end

x = cell(numel(kk), 1);
y = cell(numel(kk), 1);
z = cell(numel(kk), 1);

i = 1;
for k = kk
    x{i} = f.domain.x{k};
    y{i} = f.domain.y{k};
    z{i} = f.domain.z{k};
    i = i+1;
end

if ( nargin > 1 && numel(kk) == 1 )
    x = x{1};
    y = y{1};
    z = z{1};
end

end
