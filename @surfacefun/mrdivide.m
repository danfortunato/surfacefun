function X = mrdivide(B, A)
%\   Right matrix divide for SURFACEFUN objects.
%   B/A gives the least squares solution to X*A = B, where A, B, and X may
%   be numeric matrices or SURFACEFUN quasimatrices. For a scalar A and
%   SURFACEFUN B, B/A divides B by A.
%
%   See also MLDIVIDE, RDIVIDE.

if ( isnumeric(A) && isscalar(A) )
    X = rdivide(B, A);
elseif ( size(A, 2) ~= size(B, 2) )
    error('SURFACEFUN:mrdivide:dims', 'Matrix dimensions must agree.')
elseif ( isnumeric(A) && isa(B, 'surfacefun') )
    % [INF x M] * [M x N] = [INF x N]:
    [Q, R] = qr(B);
    X = Q * (R/A);
elseif ( isa(A, 'surfacefun') && isnumeric(B) )
    % [M x INF] * [INF x N] = [M x N]:
    [Q, R] = qr(A);
    X = (B/R) * Q';
elseif ( isa(A, 'surfacefun') && isa(B, 'surfacefun') )
    % [M x N] * [N x INF] = [M x INF]:
    [Q, R] = qr(A');
    X = (B*Q) / R';
else
    error('SURFACEFUN:mrdivide:mrdivide', 'Not supported. Did you mean ./ ?');
end

end
