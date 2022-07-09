function dom = refine(dom, nref)

if ( isempty(dom) )
    return
end

if ( nargin < 2 )
    nref = 1;
end

n = size(dom.x{1}, 1);
x  = chebpts(n, [-1 1]);
xL = chebpts(n, [-1 0]);
xR = chebpts(n, [ 0 1]);
BL = barymat(xL, x);
BR = barymat(xR, x);

x = dom.x;
y = dom.y;
z = dom.z;

for r = 1:nref
    xnew = cell(4*length(x), 1);
    ynew = cell(4*length(y), 1);
    znew = cell(4*length(z), 1);
    l = 0;
    for k = 1:length(x)
        xnew{l+1} = BL * x{k} * BL.'; % Lower left
        ynew{l+1} = BL * y{k} * BL.';
        znew{l+1} = BL * z{k} * BL.';
        xnew{l+2} = BL * x{k} * BR.'; % Lower right
        ynew{l+2} = BL * y{k} * BR.';
        znew{l+2} = BL * z{k} * BR.';
        xnew{l+3} = BR * x{k} * BL.'; % Upper left
        ynew{l+3} = BR * y{k} * BL.';
        znew{l+3} = BR * z{k} * BL.';
        xnew{l+4} = BR * x{k} * BR.'; % Upper right
        ynew{l+4} = BR * y{k} * BR.';
        znew{l+4} = BR * z{k} * BR.';
        l = l+4;
    end
    x = xnew;
    y = ynew;
    z = znew;
end

dom = surfacemesh(x, y, z);

end
