function vals = chebvals2plotvals(vals, nplotpts)
%CHEBVALS2PLOTVALS   Convert 2D Chebyshev coefficients to values at equispaced points.

persistent Eval nstored nplotptsstored

if ( nargin < 2 )
    nplotpts = 100;
end

n = size(vals, 1);
if ( isempty(Eval) || n ~= nstored || nplotpts ~= nplotptsstored )
    nstored = n;
    nplotptsstored = nplotpts;
    x = linspace(-1, 1, nplotpts).';
    xcheb = chebpts(n);
    Eval = barymat(x, xcheb);
end

vals = Eval * vals * Eval.';

end
