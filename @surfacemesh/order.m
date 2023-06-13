function p = order(dom)
%ORDER   Order of a SURFACEMESH.

if ( isempty(dom) )
    p = [];
else
    p = size(dom.x{1}, 1) - 1;
end

end
