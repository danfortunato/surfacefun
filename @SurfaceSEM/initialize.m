function initialize(S, varargin)
%INITIALIZE   Initialize and solve local subproblems.
%   INITIALIZE(S, RHS) will initialize the SURFACESEM object S on each of
%   the subpatches of S.domain using the righthand side RHS.
%
%   INITIALIZE(S) assumes the problem is homogeneous (i.e., RHS = 0).
%
%   The full sequence for solving a problem using a SURFACESEM object S
%   is:
%
%      initialize(S, RHS)
%      build(S)   % (optional)
%      sol = S\bc % or sol = solve(S, bc)
%
%   See also BUILD, SOLVE.

% Initialize all leaf patches:
S.patches = SurfaceSEM.Leaf.initialize(S.domain, varargin{:});

end
