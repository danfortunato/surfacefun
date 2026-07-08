function vals = koornvals2plotvals(vals, nplotpts)
%KOORNVALS2PLOTVALS   Convert Koornwinder coefficients to values at equispaced points.

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
    [xplot, yplot] = trianglepts(nplotpts, type='linspace');
    K = koornwinder(n-1, x, y);
    Kplot = koornwinder(n-1, xplot, yplot);
    Eval = Kplot / K;
end

vals = Eval * vals;

end
