function dom = square(n, nref, rect)
%SQUARE   Create a square Cartesian mesh.

if ( nargin < 2 )
    nref = 0;
end

if ( nargin < 3 )
    rect = [-1 1 -1 1];
end

v = rect;
for l = 1:nref
    nv = size(v, 1);
    vnew = zeros(4*nv, 4);
    for k = 1:nv
        vk = v(k,:);
        mid = [mean(vk(1:2)) mean(vk(3:4))];
        vnew((k-1)*4+(1:4),:) = [ vk(1)  mid(1) vk(3)  mid(2) ;
                                  mid(1) vk(2)  vk(3)  mid(2) ;
                                  vk(1)  mid(1) mid(2) vk(4)  ;
                                  mid(1) vk(2)  mid(2) vk(4)  ];
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

x = cell(nv, 1);
y = cell(nv, 1);
z = cell(nv, 1);
for k = 1:nv
    x{k} = uu(:,:,k);
    y{k} = vv(:,:,k);
    z{k} = 0*uu(:,:,k);
end

dom = surfacemesh(x, y, z);

end
