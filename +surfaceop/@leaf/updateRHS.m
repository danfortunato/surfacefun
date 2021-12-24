function P = updateRHS(P, rhs)
%UPDATERHS   Update RHS of an SURFACEOP.LEAF object.
%   P = UPDATERHS(P, F) replaces the existing RHS of an initialized
%   SURFACEOP.LEAF object P with that given in F, which must be a matrix
%   (or cell array) of tensor-product Chebyshev values. Since only one
%   subproblem needs to be solved on each patch in the update process
%   (rather than O(n) in the original initialization) this can lead to
%   considerable performance gains when solving for multiple RHSs.
%
%   See also SURFACEOP.LEAF.INITIALIZE.

% Developer note: At user levels this is typically called with 
%  >> P.rhs = F

n = P.n;
dom = P.domain;
id = P.id;

ii = false(n);
ii(2:n-1,2:n-1) = true;
ii([1,end],[1,end]) = true;

if ( iscell(rhs) )
    rhs = rhs{1};
end

% Define scalar RHSs:
if ( isnumeric(rhs) && isscalar(rhs) )
    % Constant RHS.
    rhs = repmat(rhs, n^2, 1);
elseif ( isnumeric(rhs) && ~isscalar(rhs) )
    % We already have the coefficients of the RHS.
elseif ( isa(rhs, 'function_handle') )
    rhs = feval(rhs, dom.x{id}, dom.y{id}, dom.z{id});
end

% Restrict to interior nodes.
rhs = rhs(ii);

% Apply A^-1
if ( isempty(P.Ainv) )
    error('SURFACEOP:LEAF:updateRHS:operatorNotStored', ...
        'Discretized operator A was not stored. Cannot update RHS.');
end
S = P.Ainv(rhs);

% Replace solution operator for corners with interp conditions:
S([1:2,end-1:end],:) = 0;

% Append boundary points to solution operator:
tmpS = zeros(n^2, 1);
tmpS(ii) = S;
S = tmpS;

P.S(:,end) = S;

% Normal derivative:
P.D2N(:,end) = P.normal_d * S;

end
