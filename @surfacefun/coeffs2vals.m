function V = coeffs2vals(C)
%COEFFS2VALS   Convert a cell array of 2D Chebyshev coefficients to values.

V = cell(size(C));
for k = 1:length(C)
    V{k} = chebtech2.coeffs2vals(chebtech2.coeffs2vals(C{k}).').';
end

end
