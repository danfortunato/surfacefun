function V = volume(dom)
%VOLUME   Compute the volume of a surface.
%
%   See also SURFACEAREA.

if ( isempty(dom) )
    V = 0;
    return
end

[nv, nu] = size(dom.x{1});
wu = chebtech2.quadwts(nu); wu = wu(:);
wv = chebtech2.quadwts(nv); wv = wv(:);

V = 0;
for k = 1:length(dom)
    I = (dom.x{k} .* dom.facenormals{k}(:,:,1) + ...
         dom.y{k} .* dom.facenormals{k}(:,:,2) + ...
         dom.z{k} .* dom.facenormals{k}(:,:,3)) / 3;
    V = V + sum(sum(I .* wv .* wu.' .* sqrt(dom.J{k})));
end

%%% Alternative way:
% x = surfacefun(dom.x, dom);
% y = surfacefun(dom.y, dom);
% z = surfacefun(dom.z, dom);
% slice = @(c,k) cellfun(@(x) x(:,:,k), c, 'UniformOutput', false);
% nx = surfacefun(slice(dom.facenormals, 1), dom);
% ny = surfacefun(slice(dom.facenormals, 2), dom);
% nz = surfacefun(slice(dom.facenormals, 3), dom);
% V = integral2((x.*nx+y.*ny+z.*nz)/3);

end
