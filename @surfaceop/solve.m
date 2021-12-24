function sol = solve(S, bc)
%SOLVE   Perform a surface solve.
%   SOL = SOLVE(S, BC) returns a cell array SOL of function values 
%   representing the solution to the PDE, subject to the boundary
%   conditions specified by the function handle BC, contained in the
%   SURFACEOP object S. If S has not yet been built (see BUILD()) then
%   SOLVE() will build it. If S has not yet been initialized (see
%   INITIALIZE()) then an error is thrown.
%
%   The full sequence for solving a problem using a SURFACEOP object S is:
%
%      initialize(S, OP, RHS)
%      build(S)
%      sol = S\bc % or sol = solve(S, bc)
%
%   See also BUILD, INITIALIZE.

% Build if required:
if ( ~isBuilt(S) )
    build(S);
end

if ( ~isInitialized(S) )
    error('The object has not been initialized.');
end

if ( nargin == 1 )
    bc = [];
end

% Solve the patch object:
sol = solve(S.patches{1}, bc);

% Package into a surfacefun:
sol = surfacefun(sol, S.domain);

end
