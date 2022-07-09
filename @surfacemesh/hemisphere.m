function dom = hemisphere(n, nref, type)

switch lower(type)
    case 'l'
        comp = @le;
    case 'r'
        comp = @ge;
    otherwise
        error('Unknown hemisphere type.');
end

dom = surfacemesh.sphere(n, nref);
i = [];
for k = 1:length(dom)
    if ( all(comp(dom.x{k}, 0), 'all') )
        i = [i k];
    end
end

dom = surfacemesh(dom.x(i), dom.y(i), dom.z(i));

end
