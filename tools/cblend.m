function map = cblend(colors, n)
%CBLEND   Make a blended colormap.
%   CBLEND(COLORS) creates a colormap by linearly blending the colors in
%   the colormap matrix COLORS. The matrix columns are interpreted as RGB
%   values in the range [0 1] and the matrix must have exactly 3 columns.
%   Blending is done for each pair of colors given in COLORS.
%
%   CBLEND(COLORS, N) uses N colors to blend between each pair of colors.

%   Copyright 2021 Dan Fortunato.

if ( nargin == 1 )
    n = 100;
end

m = size(colors, 1);
map = zeros((n+1)*(m-1)+1, 3);
for k = 1:m-1
    c1 = colors(k,:);
    c2 = colors(k+1,:);
    idx = (n+1)*(k-1) + (1:n+2);
    map(idx, 1) = linspace(c1(1), c2(1), n+2);
    map(idx, 2) = linspace(c1(2), c2(2), n+2);
    map(idx, 3) = linspace(c1(3), c2(3), n+2);
end

end
