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
I = zeros(length(f), size(f, 2));
sz1 = cellfun(@(x) size(x,1), f(:,1).vals);
sz2 = cellfun(@(x) size(x,2), f(:,1).vals);
if ( all(sz1 == sz1(1)) && all(sz2 == sz2(1)) )
    % All elements are the same order, so they can share quadrature weights.
    % QUADWTS handles both triangle and quad patches and already includes
    % the square root of the Jacobian on each element.
    W = quadwts(f(:,1).domain);
    for k = 1:size(f, 2)
        for j = 1:length(f(:,k))
            I(j,k) = sum(f(:,k).vals{j} .* W(:,:,j), 'all');
        end
    end
else
    % Elements need different quadrature weights.
    for k = 1:size(f, 2)
        for j = 1:length(f(:,k))
            [nv, nu] = size(f(:,k).vals{j});
            wu = chebtech2.quadwts(nu); wu = wu(:);
            wv = chebtech2.quadwts(nv); wv = wv(:);
            I(j,k) = sum(sum(f(:,k).vals{j} .* wv .* wu.' .* sqrt(f(:,k).domain.J{j})));
        end
    end
end

% Combine norms on each element.
if ( reduce )
    I = sum(I);
end

end
