function W = quadwts(dom)
%QUADWTS   Quadrature weights for a SURFACEMESH.
%   W = QUADWTS(DOM) computes quadrature weights for the SURFACEMESH DOM.
%   These weights are computed as the native tensor-product Clenshaw-
%   Curtis quadrature weights scaled by the square root of the Jacobian on
%   each element. If DOM uses an NV x NU discretization on each element,
%   then W is a tensor of size NV x NU x LENGTH(DOM).

if ( isempty(dom) )
    W = [];
    return
end

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N points.
nelem = length(dom);
[nv, nu] = size(dom.x{1});
wu = chebtech2.quadwts(nu); wu = wu(:);
wv = chebtech2.quadwts(nv); wv = wv(:);
W = zeros(nv, nu, nelem);
for k = 1:nelem
    W(:,:,k) = wv .* wu.' .* sqrt(dom.J{k});
end

end
