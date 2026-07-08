function domnew = quad2tri(dom)
%QUAD2TRI   Convert a quadrilateral mesh to a triangle mesh.

triIdx  = ( dom.ptype == surfacemesh.patchtype.tri  );
quadIdx = ( dom.ptype == surfacemesh.patchtype.quad );

ntri  = sum(triIdx);
nquad = sum(quadIdx);

xnew = cell(ntri + 2*nquad, 1);
ynew = cell(ntri + 2*nquad, 1);
znew = cell(ntri + 2*nquad, 1);

x = dom.x;
y = dom.y;
z = dom.z;

if ( quadIdx(1) )
    n = size(x{1}, 1);
else
    npts = length(x{1});
    n = (sqrt(8*npts+1)-1) / 2;
end

[xtri, ytri] = trianglepts(n, domain='biunit');
V2C = chebtech2.vals2coeffs(eye(n));
V2C = kron(V2C, V2C);
T1 = chebpoly2(n, xtri, ytri);
B1 = T1 * V2C;
T2 = chebpoly2(n, -xtri, -ytri);
B2 = T2 * V2C;

knew = 1;
for k = 1:length(dom)
    if ( triIdx(k) )
        xnew(knew) = x(k);
        ynew(knew) = y(k);
        znew(knew) = z(k);
        knew = knew+1;
    else
        xnew{knew}   = B1*x{k}(:);
        ynew{knew}   = B1*y{k}(:);
        znew{knew}   = B1*z{k}(:);
        xnew{knew+1} = B2*x{k}(:);
        ynew{knew+1} = B2*y{k}(:);
        znew{knew+1} = B2*z{k}(:);
        knew = knew+2;
    end
end

domnew = surfacemesh(xnew, ynew, znew, 'tri');

end

function V = chebpoly2(n, x, y)

x = x(:);
y = y(:);
basis1d = chebpoly(0:n-1);
xbasis = basis1d(x);
ybasis = basis1d(y);
V = zeros(length(x), n^2);
for i = 1:n
    for j = 1:n
        V(:,(i-1)*n+j) = xbasis(:,i) .* ybasis(:,j);
    end
end

end

function vals = bary2d(chebvals, x, y)

vals = 0*x;
C = chebtech2.vals2coeffs( chebtech2.vals2coeffs( chebvals ).' ).';
Cy = chebtech2.clenshaw(y(:), C).';
for k = 1:numel(x)
    vals(k) = chebtech2.clenshaw(x(k), Cy(:,k));
end

end
