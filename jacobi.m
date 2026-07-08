function p = jacobi(n, a, b, x)
%JACOBI   Evaluate Jacobi polynomials.
%   P = JACOBI(N, A, B, X) evaluates the Jacobi polynomials P_N^{A,B} at the
%   points X, for N >= 0, A, B > -1, and X in [-1, 1].

n = n(:);
x = x(:);

if ( any(n < 0) )
    error('N must be greater than 0.');
end

if ( a <= -1 || b <= -1 )
    error('A and B must be greater than -1.');
end

nmax = max(n);
m = numel(x);
p = zeros(m, nmax+1);

p(:,1) = 1;
if ( nmax < 1 ), return, end
p(:,2) = (1+(a+b)/2)*x  + (a-b)/2;
for k = 2:nmax
    c1 = 2*k*(k+a+b)*(2*k-2+a+b);
    c2 = (2*k-1+a+b)*(2*k+a+b)*(2*k-2+a+b);
    c3 = (2*k-1+a+b)*(a+b)*(a-b);
    c4 = -2*(k-1+a)*(k-1+b)*(2*k+a+b);
    p(:,k+1) = ((c3+c2*x).*p(:,k) + c4*p(:,k-1)) / c1;
end
p = p(:,n+1);

end
