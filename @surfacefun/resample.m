function u = resample(u, n)

if ( isempty(u) )
    return
end

m = size(u.vals{1}, 1);
B = barymat(chebpts(n), chebpts(m));

vals = u.vals;
for k = 1:length(u)
    vals{k} = B * vals{k} * B.';
end

dom = resample(u.domain, n);
u = surfacefun(vals, dom);

end
