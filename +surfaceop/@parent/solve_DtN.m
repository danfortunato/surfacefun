function u = solve_DtN(P, bc)
%SOLVE   Solve a parent patch.
%   U = SOLVE(P, BC) returns a cell array U of function values representing
%   the PDE solution on the parent P with Dirichlet boundary data given by
%   BC.

if ( ~isnumeric(bc) )
    % Evaluate the RHS if given a function handle:
    bc = feval(bc, P.xyz(:,1), P.xyz(:,2), P.xyz(:,3));
    bc = repmat(bc, 1, size(P.u_part, 2));
elseif ( isscalar(bc) )
    % Convert a scalar to a constant vector:
    bc = repmat(bc, size(P.xyz, 1), size(P.u_part, 2));
end

% Evaluate the solution operator for the parent:
u = P.u_part;
if ( ~isempty(bc) )
    u = u + P.S * bc;
end

% Construct boundary conditions for children and recurse.

% Construct boundary indices:
i1 = 1:numel(P.idx1{1});
i2 = 1:numel(P.idx2{1});
if ( ~isempty(i1) )
    i2 = i2 + i1(end);
end
% CAT() is 10x faster than CELL2MAT().
idx1 = cat(1, P.idx1{:}); % idx1 = cell2mat(P.idx1.');
idx2 = cat(1, P.idx2{:}); % idx2 = cell2mat(P.idx2.');

% Assemble boundary conditions for child patches:
ubc1 = ones(size(P.child1.S, 2)-1, size(P.u_part, 2));
ubc1(idx1,:) = [bc(i1,:) ; P.flip1.'*u];
ubc2 = ones(size(P.child2.S, 2)-1, size(P.u_part, 2));
ubc2(idx2,:) = [bc(i2,:) ; P.flip2.'*u];

% Solve for the child patches:
u1 = solve_DtN(P.child1, ubc1);
u2 = solve_DtN(P.child2, ubc2);

% Concatenate for output:
u = [u1 ; u2]; 

end
