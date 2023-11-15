function X = mldivide(A, B)
%\   Left matrix divide for SURFACEFUN objects.
%   A\B gives the least squares solution to A*X = B, where A, B, and X may
%   be numeric matrices or SURFACEFUN quasimatrices. For a scalar A and
%   SURFACEFUN B, A\B divides B by A.
%
%   See also MRDIVIDE, LDIVIDE.

if ( isnumeric(A) && isscalar(A) )
    X = ldivide(A, B);
elseif ( size(A, 1) ~= size(B, 1) )
    error('SURFACEFUN:mldivide:dims', 'Matrix dimensions must agree.')
elseif ( isnumeric(A) && isa(B, 'surfacefun') )
    % [M x N] * [N x INF] = [M x INF]:
    [Q, R] = qr(B');
    X = (A\R') * Q';
elseif ( isa(A, 'surfacefun') && isnumeric(B) )
    % [M x INF] * [INF x N] = [M x N]:
    [Q, R] = qr(A');
    X = Q * (R'\B);
elseif ( isa(A, 'surfacefun') && isa(B, 'surfacefun') )
    % [INF x N] * [N x M] = [INF x M]:
    [Q, R] = qr(A);
    X = R \ (Q'*B);
else
    error('SURFACEFUN:mldivide:mldivide', 'Not supported. Did you mean .\ ?');
end

end
