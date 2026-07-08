function vals = koornvals2duffyvals(vals, nplotpts)
%KOORNVALS2PLOTVALS   Convert Koornwinder coefficients to values at equispaced
% tensor-product points.

persistent Eval nstored nplotptsstored

if ( nargin < 2 )
    nplotpts = 100;
end

npts = length(vals);
n = (sqrt(8*npts+1)-1) / 2;
if ( isempty(Eval) || n ~= nstored || nplotpts ~= nplotptsstored )
    nstored = n;
    nplotptsstored = nplotpts;
    [x, y] = trianglepts(n);
    [eta1, eta2] = meshgrid(linspace(0, 1, nplotpts));
    % Duffy transformation
    xplot = eta1.*(1-eta2);
    yplot = eta2;
    K = koornwinder(n-1, x, y);
    Kplot = koornwinder(n-1, xplot, yplot);
    Eval = Kplot / K;
end

vals = Eval * vals;

end
