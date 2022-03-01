function dom = resample(dom, n)

if ( isempty(dom) )
    return
end

m = size(dom.x{1}, 1);
B = barymat(chebpts(n), chebpts(m));

x = dom.x;
y = dom.y;
z = dom.z;
for k = 1:length(dom)
    x{k} = B * x{k} * B.';
    y{k} = B * y{k} * B.';
    z{k} = B * z{k} * B.';
end

dom = surfacemesh(x, y, z);

end
