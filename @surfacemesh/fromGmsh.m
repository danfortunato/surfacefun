function dom = fromGmsh(msh)

names = fields(msh);
field = names(startsWith(names, 'QUADS'));
field = field{1};
quads = msh.(field);

nelem = size(quads, 1);
n = sqrt(size(quads, 2) - 1);
x = cell(nelem, 1);
y = cell(nelem, 1);
z = cell(nelem, 1);

for k = 1:nelem
    x{k} = quadTensorProductOrder(msh.POS(quads(k, 1:n^2), 1), n);
    y{k} = quadTensorProductOrder(msh.POS(quads(k, 1:n^2), 2), n);
    z{k} = quadTensorProductOrder(msh.POS(quads(k, 1:n^2), 3), n);
end

ex = linspace(-1, 1, n).';
cheb2equi = bary(ex, eye(n));
equi2cheb = inv(cheb2equi);
for k = 1:nelem
    x{k} = equi2cheb * x{k} * equi2cheb.';
    y{k} = equi2cheb * y{k} * equi2cheb.';
    z{k} = equi2cheb * z{k} * equi2cheb.';
end

dom = surfacemesh(x, y, z);

end

function y = quadTensorProductOrder(x, n)

y = zeros(n);

y(1) = x(1);

if ( n > 1 )
    y(n) = x(2);
    y(n^2) = x(3);
    y(n*(n-1)+1) = x(4);
end

for i = 1:n-2
    y(i+1) = x(4+i);
end

for i = 1:n-2
    y(n^2-i) = x(4+2*(n-2)+i);
end

% Interior (recursion)
if ( n > 2 )
    % Reorder interior nodes
    interior = x((4+4*(n-2)+1):end);
    interior = quadTensorProductOrder(interior, n-2);

    % Copy into full node list
    for j = 1:n-2
        y(j*n+1) = x(4+3*(n-2)+n-2-j+1);
        for i = 1:n-2
            y(j*n+i+1) = interior((j-1)*(n-2)+(i-1)+1);
        end
        y((j+1)*n) = x(4+(n-2)+j);
    end
end

end
