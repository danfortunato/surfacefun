function [Q, R] = qr(f)
%QR   QR factorization of an array-valued SURFACEFUN.
%   [Q, R] = QR(F) returns a QR factorization of F such that F = Q*R, where
%   the array of SURFACEFUNs Q is orthogonal (with respect to the
%   continuous L^2 norm on the surface) and of the same size as F and R is
%   an M x M upper-triangular matrix when F has M columns.

% If f has only one column we simply scale it:
if ( size(f, 2) == 1 )
    R = sqrt(integral2(f.^2));
    Q = f ./ R;
    return
end

dom = f(:,1).domain;
A = zeros(numel(dom), size(f, 2));
for k = 1:size(f, 2)
    A(:,k) = reshape(cell2mat(f(:,k).vals.'), [], 1);
end

I = quadwts(dom);
M = spdiags(I(:), 0, numel(dom), numel(dom));
W = sqrt(M);       % In general, chol(M)
[U, R] = qr(W*A, 0);
s = sign(diag(R)); % }
s(~s) = 1;         % } Enforce diag(R) >= 0
S = diag(s);       % }
QQ = (W \ U) * S;
R = S * R;

% Convert to surfacefuns:
Q = f;
nelem = length(dom);
[nv, nu] = size(dom.x{1});
for k = 1:size(f, 2)
    vals = reshape(QQ(:,k), nv, nu, nelem);
    vals = squeeze(mat2cell(vals, nv, nu, ones(nelem, 1)));
    Q(:,k) = surfacefun(vals, dom);
end

end
