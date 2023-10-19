classdef leaf < surfaceop.patch
%SURFACEOP.LEAF   Leaf subclass of a patch (where subproblems are solved).
%   P = SURFACEOP.LEAF(DOMAIN, S, D2N, EDGES, AINV) creates a LEAF object
%   P and assigns each of the inputs to their associated properties in P.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTES:
% * Note: We assume that a LEAF has four sides and that its boundary nodes
%   XYZ are stored in the order "left", "right", "down", "up".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties

        n        % Discretization size of patch.
        Ainv     % Local (homogeneous BC) solution operator.
        normal_d % Normal derivative operator.
                 % (These are stored so the RHS can be efficiently updated)

    end

    methods

        function P = leaf(dom, n, id, S, BtB, BtB_scl, u_part, du_part, edges, xyz, w, Ainv, normal_d)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = dom;           % Domain.
            P.n = n;                  % Discretization size.
            P.id = id;                % Index of patch in domain.
            P.S = S;                  % Solution operator.
            P.BtB = BtB;              % Poincare-Steklov operator.
            P.BtB_scl = BtB_scl;      % Scalings for Poincare-Steklov operator.
            P.u_part = u_part;        % Particular solution.
            P.du_part = du_part;      % Outgoing boundary data from particular solution.
            P.edges = edges;          % Boundary edges.
            P.xyz = xyz;              % Boundary nodes.
            P.w = w;                  % Boundary quadrature weights.
            P.Ainv = Ainv;            % Local solution operator.
            P.normal_d = normal_d;    % Normal derivative operator.
            P.len = 1;

        end

    end

    methods ( Static )

        % Initialize an array of LEAF objects.
        P = initialize_DtN(op, dom, rhs);
        P = initialize_ItI(op, dom, eta, rhs);

    end

end
