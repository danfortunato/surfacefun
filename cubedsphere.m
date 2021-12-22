function [x, y, z] = cubedsphere(n, nref, type)

if ( nargin == 0 )
    test_cubedsphere();
    return
end

if ( nargin < 2 )
    nref = 0;
end

if ( nargin < 3 )
    type = 'quasiuniform';
end

switch lower(type)
    case 'naive'
        project = @projectNaive;
    case 'quasiuniform'
        project = @projectQuasiUniform;
    otherwise
        error('Unknown projection type.');
end

v = [-1 1 -1 1];
for l = 1:nref
    nv = size(v, 1);
    vnew = zeros(4*nv, 4);
    for k = 1:nv
        vk = v(k,:); 
        mid = [mean(vk(1:2)) mean(vk(3:4))];
        vnew((k-1)*4+(1:4),:) = [ vk(1)  mid(1) vk(3)  mid(2) ;
                                  mid(1) vk(2)  vk(3)  mid(2) ;
                                  mid(1) vk(2)  mid(2) vk(4)  ;
                                  vk(1)  mid(1) mid(2) vk(4)  ];
    end
    v = vnew;
end

[xx0, yy0] = chebpts2(n, n, [0 1 0 1]);
nv = size(v, 1);
uu = zeros(n, n, nv);
vv = zeros(n, n, nv);
for k = 1:nv
    sclx = diff(v(k,1:2));
    scly = diff(v(k,3:4));
    uu(:,:,k) = sclx*xx0 + v(k,1);
    vv(:,:,k) = scly*yy0 + v(k,3);
end

xx = cat(3, -1+0*uu, 1+0*uu, vv, vv, uu, uu);
yy = cat(3, uu, uu, -1+0*uu, 1+0*uu, vv, vv);
zz = cat(3, vv, vv, uu, uu, -1+0*uu, 1+0*uu);
[xx, yy, zz] = project(xx, yy, zz);

x = cell(6*nv, 1);
y = cell(6*nv, 1);
z = cell(6*nv, 1);
for k = 1:6*nv
    x{k} = xx(:,:,k);
    y{k} = yy(:,:,k);
    z{k} = zz(:,:,k);
end

end

function [xp, yp, zp] = projectNaive(x, y, z)
    nrm = sqrt(x.^2 + y.^2 + z.^2);
    xp = x ./ nrm;
    yp = y ./ nrm;
    zp = z ./ nrm;
end

function [xp, yp, zp] = projectQuasiUniform(x, y, z)
    xp = x .* sqrt(1 - y.^2/2 - z.^2/2 + (y.^2.*z.^2)/3);
    yp = y .* sqrt(1 - x.^2/2 - z.^2/2 + (x.^2.*z.^2)/3);
    zp = z .* sqrt(1 - x.^2/2 - y.^2/2 + (x.^2.*y.^2)/3);
end

function test_cubedsphere()

clf

% Plot a white sphere
[x, y, z] = sphere(100);
scl = 0.99;
surf(scl*x, scl*y, scl*z, 1+0*x, 'FaceColor', 'w', 'EdgeColor', 'none')
axis equal
hold on

n = 20;
nref = 3;
[x, y, z] = cubedsphere(n, nref, 'quasiuniform');
ne = length(x);
x = cat(3, x{:});
y = cat(3, y{:});
z = cat(3, z{:});
X = [ x(:,1,:) ; nan(1,1,ne) ; x(:,n,:) ; nan(1,1,ne) ; permute(x(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(x(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
Y = [ y(:,1,:) ; nan(1,1,ne) ; y(:,n,:) ; nan(1,1,ne) ; permute(y(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(y(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
Z = [ z(:,1,:) ; nan(1,1,ne) ; z(:,n,:) ; nan(1,1,ne) ; permute(z(1,:,:), [2 1 3]) ; nan(1,1,ne) ; permute(z(n,:,:), [2 1 3]) ; nan(1,1,ne) ];
plot3(X(:), Y(:), Z(:), 'k-', 'LineWidth', 1)
shg

end

