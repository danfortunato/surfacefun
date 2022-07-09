function u = refine(u, nref)

if ( isempty(u) )
    return
end

if ( nargin < 2 )
    nref = 1;
end

dom = u.domain;
n = size(dom.x{1}, 1);
x  = chebpts(n, [-1 1]);
xL = chebpts(n, [-1 0]);
xR = chebpts(n, [ 0 1]);
BL = barymat(xL, x);
BR = barymat(xR, x);

vals = u.vals;
for r = 1:nref
    rvals = cell(4*length(vals), 1);
    l = 0;
    for k = 1:length(vals)
        rvals{l+1} = BL * vals{k} * BL.'; % Lower left
        rvals{l+2} = BL * vals{k} * BR.'; % Lower right
        rvals{l+3} = BR * vals{k} * BL.'; % Upper left
        rvals{l+4} = BR * vals{k} * BR.'; % Upper right
        l = l+4;
    end
    vals = rvals;
end

rdom = refine(dom, nref);
u = surfacefun(vals, rdom);

end
