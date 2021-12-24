function [x, y] = plotpts(f, kk)
%PLOTPTS   Get plot points.
%   [X, Y] = PLOTPTS(F) returns points (X, Y) which are uniformly spaced
%   over each patch in SOL, to be used for plotting.
%
%   [X, Y] = PLOTPTS(F, KK) returns plot points for the KK-th patches.

d = f.domain;
vals = f.vals;

x = cell(size(vals));
y = cell(size(vals));

if ( nargin < 2 )
    kk = 1:size(d, 1);
end

kk = kk(:).';

for k = kk
    if ( isnumeric(d(k,:)) )
        [x{k,1}, y{k,1}] = meshgrid(linspace(d(k,1), d(k,2), f.nplotpts), ...
                                    linspace(d(k,3), d(k,4), f.nplotpts));
    else
        [xk, yk] = meshgrid(linspace(-1, 1, f.nplotpts));
        [x{k,1}, y{k,1}] = transformGrid(d(k,:), xk, yk);
    end
end

if ( nargin > 1 && numel(kk) == 1 )
    x = x{1};
    y = y{1};
end

end
