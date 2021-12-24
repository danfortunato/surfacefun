function vals = chebvals2plotvals(vals)
%CHEBVALS2PLOTVALS   Convert 2D Chebyshev coefficients to values at equispaced points.

persistent Eval nstored
nplotpts = 100;

n = size(vals, 1);
if ( isempty(Eval) || n ~= nstored )
    nstored = n;
    x = linspace(-1, 1, nplotpts).';
    xcheb = chebpts(n);
    Eval = barymat(x, xcheb);
end

vals = Eval * vals * Eval.';

end
