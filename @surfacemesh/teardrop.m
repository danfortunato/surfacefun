function dom = teardrop(n, nref, type)
%TEARDROP   Create a teardrop mesh.

if ( nargin < 2 )
    nref = 0;
end

if ( nargin < 3 )
    type = 'naive';
end

switch lower(type)
    case 'naive'
        project = @projectNaive;
    case 'quasiuniform'
        project = @projectQuasiUniform;
    otherwise
        error('SURFACEMESH:TEARDROP:projection', 'Unknown projection type.');
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

xx = cat(3, -1+0*uu, 1+0*uu, vv,      uu,     uu,      vv);
yy = cat(3, uu,      vv,     -1+0*uu, 1+0*uu, vv,      uu);
zz = cat(3, vv,      uu,     uu,      vv,     -1+0*uu, 1+0*uu);
[xx, yy, zz] = project(xx, yy, zz);
[xx, yy, zz] = project_teardrop(xx, yy, zz);

x = cell(6*nv, 1);
y = cell(6*nv, 1);
z = cell(6*nv, 1);
for k = 1:6*nv
    x{k} = xx(:,:,k);
    y{k} = yy(:,:,k);
    z{k} = zz(:,:,k);
end

dom = surfacemesh(x, y, z);

end

function [xp, yp, zp] = project_teardrop(x, y, z)

m = 1;
th = cart2pol(x, y);
t = acos(z);
rnew = sin(t).*sin(t/2).^m;
[xp, yp, zp] = pol2cart(th, rnew, z);

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
