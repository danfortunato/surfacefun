classdef ( Abstract ) patch
%SURFACEOP.PATCH   Abstract patch object.

    properties

        domain  % Domain of patch.
        id      % Index of patch in domain.
        S       % Solution operator for patch.
        BtB     % Poincare-Steklov operator for patch.
        BtB_scl % Cell array of scalars or function handles.
                % The k-th entry is the scaling for the BtB map on side k.
        u_part  % Particular solution.
        du_part % Outgoing boundary data from particular solution.
        edges   % Boundary edges of patch.
        xyz     % Boundary grid points of patch.
        w       % Boundary quadrature weights of patch.
        len

    end

    methods ( Abstract )

        % Solve a patch.
        [u, d] = solve_DtN(P, bc);
        [u, d] = solve_ItI(P, bc);

        % Update RHS of a patch.
        f = updateRHS(f, rhs);

        % Number of degrees of freedom in a patch.
        N = numel(P);

        % Number of patches in a patch.
        N = length(P);

    end

end
