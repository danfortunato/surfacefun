function A = morton(m, n)
%MORTON   Create a Morton ordering matrix.

if ( nargin == 1 )
    n = m;
end

if ( m == 0 || n == 0 )
    A = [];
    return
end

isPowerM = bitand(m, m-1) == 0;
isPowerN = bitand(n, n-1) == 0;
if ( ~isPowerM || ~isPowerN )
    error('M and N must be a power of 2.');
end

if ( m == 1 )
    A = 1:n;
elseif ( n == 1 )
    A = (1:m).';
elseif ( m == 2 && n == 2 )
    A = [1 2; 3 4];
else
    B = morton(m/2, n/2);
    A = [B B+(m/2)*(n/2); B+(m/2)*(n/2)*2 B+(m/2)*(n/2)*3];
end

end
