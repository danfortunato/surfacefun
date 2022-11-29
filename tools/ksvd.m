function [U, S, V] = ksvd(A, m, n, r)
%KSVD   Kronecker singular value decomposition.
%   [U, S, V] = KSVD(A, M, N, R) produces cell arrays U and V of length R
%   and a vector S of length R such that
%
%      || A - sum_k S(k) kron(U{k}, V{k}) ||_F = min,
%
%   where U{k} is M x M and V{k} is N x N.
%
%   [U, S, V] = KSVD(A, M, N) uses R = MIN(M, N).
%
%   [U, S, V] = KSVD(A, N) uses M = N.

if ( nargin == 0 )
    test_ksvd();
    return
end

if ( nargin < 3 )
    n = m;
end

if ( nargin < 4 )
    r = min(m, n);
end

mn = size(A, 1);
if ( mn ~= size(A, 2) )
    error('Input must be a square matrix.');
end

if ( mn ~= m*n )
    error('The product of the block dimensions must equal the matrix size.');
end

% If we write the matrix A as an M x M block matrix with N x N blocks,
%
%    A = [ A11 . . . A1M
%           .  .      .
%           .    .    .
%           .      .  .
%          AM1 . . . AMM ],
% 
% then the Kronecker SVD can be computed by taking the SVD of the reshaped
% M^2 x N^2 matrix
%
%    AA = [ A11(:).'
%           A21(:).'
%           A31(:).'
%            ...    ].

AA = zeros(m^2, n^2);
i = 1;
for j = 1:m
    for k = 1:m
        block = A((k-1)*n + (1:n), (j-1)*n + (1:n));
        AA(i,:) = block(:).';
        i = i + 1;
    end
end

if ( issparse(A) )
    AA = sparse(AA);
end

[UU, SS, VV] = svds(AA, r);

S = zeros(r, 1);
U = cell(r, 1);
V = cell(r, 1);
for k = 1:r
    S(k) = SS(k,k);
    U{k} = reshape(UU(:,k), m, m);
    V{k} = reshape(VV(:,k), n, n);
    if ( issparse(A) )
        U{k} = sparse(U{k});
        V{k} = sparse(V{k});
    end
end

end

function test_ksvd()

% Kronecker product sizes
m = 10;
n = 20;

% Random singular values
r = min(m^2, n^2);
A = rand(m*n);
[U, S, V] = ksvd(A, m, n, r);
B = kronsum(U, S, V);
norm(A - B)

% Exponentially decaying singular values
r = 3;
S = 10.^linspace(0, -16, r).';
A = zeros(m*n);
for k = 1:r
    A = A + S(k) * kron(rand(m), rand(n));
end
[U, S, V] = ksvd(A, m, n, r);
B = kronsum(U, S, V);
norm(A - B)

end

function A = kronsum(U, S, V)

m = size(U{1}, 1);
n = size(V{1}, 1);
r = size(S, 1);

A = zeros(m*n);
for k = 1:r
    A = A + S(k) * kron(U{k}, V{k});
end

end
