function [nu, nv] = size(dom)
%SIZE   Size of a SURFACEMESH.

if ( isempty(dom) )
    nu = [];
    nv = [];
    return
end

if ( ~all(dom.ptype == dom.ptype(1)) )
    error('Heterogeneous patch types are not yet supported.');
end

switch ( dom.ptype(1) )
    case 'tri'
        nu = 1;
        nv = length(dom.x{1});
    case 'quad'
        [nv, nu] = size(dom.x{1});
end

end
