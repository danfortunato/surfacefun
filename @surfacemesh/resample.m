function dom = resample(dom, n)

if ( isempty(dom) )
    return
end

if ( ~all(dom.ptype == dom.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

x = dom.x;
y = dom.y;
z = dom.z;

switch ( dom.ptype(1) )
    case 'tri'
        m = order(dom)+1;
        [xm, ym] = trianglepts(m);
        [xn, yn] = trianglepts(n);
        Km  = koornwinder(m-1, xm, ym);
        Knm = koornwinder(m-1, xn, yn);
        B = Knm / Km;
        for k = 1:length(dom)
            x{k} = B * x{k};
            y{k} = B * y{k};
            z{k} = B * z{k};
        end

    case 'quad'
        m = size(dom.x{1}, 1);
        B = barymat(chebpts(n), chebpts(m));
        for k = 1:length(dom)
            x{k} = B * x{k} * B.';
            y{k} = B * y{k} * B.';
            z{k} = B * z{k} * B.';
        end
end

ctvy = dom.connectivity;
dom = surfacemesh(x, y, z, dom.ptype);
dom.connectivity = ctvy;

end
