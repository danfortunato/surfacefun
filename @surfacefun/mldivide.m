function X = mldivide(A, B)
%\   Left matrix divide for CHEBFUN objects.
%   For array-valued SURFACEFUNs A and B, A\B computes the least squares
%   solution X to A*X = B. For a scalar A and SURFACEFUN B, A\B divides B
%   by A.

if ( isnumeric(A) && isscalar(A) )
    X = ldivide(A, B);

elseif ( isa(A, 'surfacefun') && isa(B, 'surfacefun') )
    % (Inf x N) * (N x M) = (Inf x M):
    n = size(A, 2);
    m = size(B, 2);

    [Q, R] = qr(A);

    % Compute the L^2 inner product of each column of Q with each column
    % of B
    QB = zeros(n, m);
    for j = 1:n
        for k = 1:m
            QB(j,k) = sum2( Q(:,j) .* B(:,k) );
        end
    end

    X = R \ QB;

else
    error('SURFACEFUN:mldivide:mldivide', 'Not supported. Did you mean .\ ?');

end

end
