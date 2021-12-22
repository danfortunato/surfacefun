classdef Leaf < SurfaceSEM.Patch
%SURFACESEM.LEAF   Leaf subclass of a patch (where subproblems are solved).
%   P = SURFACESEM.LEAF(DOMAIN, S, D2N, EDGES, AINV) creates a LEAF object
%   P and assigns each of the inputs to their associated properties in P.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * Note: We assume that a LEAF has four sides and that its boundary nodes
%   XYZ are stored in the order "left", "right", "down", "up".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        n        % Discretization size of patch.
        Ainv     % Local (homogeneous BC) solution operator.
        normal_d % Normal derivative operator.
                 % (These are stored so the RHS can be efficiently updated)

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function P = Leaf(dom, S, D2N, edges, xyz, Ainv, normal_d)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.n = size(dom.x, 1);      % Discretization size.
            P.domain = dom;            % Domain.
            P.S = S;                   % Solution operator.
            P.D2N = D2N;               % Dirichlet-to-Neumann map.
            P.edges = edges;           % Edges.
            P.xyz = xyz;               % Boundary nodes.
            if ( nargin > 5 )
                P.Ainv = Ainv;         % Local solution operator.
            end
            if ( nargin > 6 )
                P.normal_d = normal_d; % Normal derivative operator.
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static )

        % Initialize an array of LEAF objects.
        P = initialize(dom, rhs);

    end

end
