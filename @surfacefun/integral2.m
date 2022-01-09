function I = integral2(f, varargin)
%INTEGRAL2   Double integral of a SURFACEFUN.
%   I = INTEGRAL2(F) returns the double integral of the SURFACEFUN F over
%   its domain.
%
%   I = INTEGRAL2(F, 'all') returns an array of double integrals over each
%   patch of F.
%
%   See also SUM2.

% Parse arguments.
reduce = true;
if ( nargin == 2 )
    if ( strcmp(varargin{1}, 'all') )
        reduce = false;
    end
end

% Empty check:
if ( isempty(f) )
    I = 0;
    return
end

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N points.
I = zeros(length(f), 1);
for k = 1:length(f)
    [nv, nu] = size(f.vals{k});
    wu = chebtech2.quadwts(nu); wu = wu(:);
    wv = chebtech2.quadwts(nv); wv = wv(:);
    I(k) = sum(sum(f.vals{k} .* wv .* wu.' .* sqrt(f.domain.J{k})));
end

% Combine norms on each element.
if ( reduce )
    I = sum(I);
end

end
