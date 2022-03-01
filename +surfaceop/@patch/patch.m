classdef ( Abstract ) patch
%SURFACEOP.PATCH   Abstract patch object.

    properties

        domain  % Domain of patch.
        id      % Index of patch in domain.
        S       % Solution operator for patch.
        D2N     % Dirichlet-to-Neumann map for patch.
        D2N_scl % Cell array of scalars or function handles.
                % The k-th entry is the scaling for the D2N map on side k.
        u_part
        du_part
        edges   % Boundary edges of patch.
        xyz     % Boundary grid points of patch.
        len

    end

    methods ( Abstract )

        % Solve a patch.
        [u, d] = solve(P, bc);

        % Update RHS of a patch.
        f = updateRHS(f, rhs);

        % Number of degrees of freedom in a patch.
        N = numel(P);

        % Number of patches in a patch.
        N = length(P);

    end

end
