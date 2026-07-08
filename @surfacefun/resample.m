function u = resample(u, n)

if ( isempty(u) )
    return
end

if ( ~all(u.domain.ptype == u.domain.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

switch ( u.domain.ptype(1) )
    case 'tri'
        m = order(u)+1;
        [xm, ym] = trianglepts(m);
        [xn, yn] = trianglepts(n);
        Km  = koornwinder(m-1, xm, ym);
        Knm = koornwinder(m-1, xn, yn);
        B = Knm / Km;
        vals = u.vals;
        for k = 1:length(u)
            vals{k} = B * vals{k};
        end

    case 'quad'
        m = size(u.vals{1}, 1);
        B = barymat(chebpts(n), chebpts(m));
        vals = u.vals;
        for k = 1:length(u)
            vals{k} = B * vals{k} * B.';
        end
end

dom = resample(u.domain, n);
u = surfacefun(vals, dom);

end
