classdef ( Abstract ) Patch
%SURFACESEM.PATCH   Abstract patch object.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain  % Domain of patch.
        S       % Solution operator for patch.
        D2N     % Dirichlet-to-Neumann map for patch.
        edges   % Boundary edges of patch.
        xyz     % Boundary grid points of patch.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Abstract )

        % Solve a patch.
        [u, d] = solve(P, bc);

        % Number of degrees of freedom in a patch.
        N = numel(P);

        % Number of patches in a patch.
        N = length(P);

    end

end
