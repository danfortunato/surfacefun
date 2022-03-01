function dom = fromRhino(filename, n)

if ( nargin < 2 )
    n = 8;
end

A = readmatrix(filename);
nelem = size(A,1)/n^2;
x = mat2cell(reshape(A(:,1), n, []), n, n*ones(nelem,1)).';
y = mat2cell(reshape(A(:,2), n, []), n, n*ones(nelem,1)).';
z = mat2cell(reshape(A(:,3), n, []), n, n*ones(nelem,1)).';
dom = surfacemesh(x, y, z);

end
