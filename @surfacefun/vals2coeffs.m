function C = vals2coeffs(V)
%VALS2COEFFS   Convert a cell array of values to 2D Chebyshev coefficients.

C = cell(size(V));
for k = 1:length(V)
    C{k} = chebtech2.vals2coeffs(chebtech2.vals2coeffs(V{k}).').';
end

end
