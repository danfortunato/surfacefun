classdef parent < surfaceop.patch
%SURFACEOP.PARENT   Parent subclass of a patch.
%   P = SURFACEOP.PARENT(DOMAIN, S, D2N, EDGES, XY, CHILD1, CHILD2, IDX1,
%   IDX2, L2G1, L2G2) creates a SURFACEOP.PARENT object P and assigns each
%   of the inputs to their associated properties in P.

    properties

        child1 = [] % Child patch
        child2 = [] % Child patch
        idx1        % How p.xyz relates to p.child1.xyz
        idx2        % How p.xyz relates to p.child2.xyz
        dA          % Decomposition of interface linear system

    end

    methods

        function P = parent(domain, id, S, D2N, dA, edges, xyz, child1, child2, idx1, idx2)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = domain;
            P.id = id;
            P.S = S;
            P.D2N = D2N;
            P.dA = dA;
            P.edges = edges;
            P.xyz = xyz;

            % Assign children:
            P.child1 = child1;
            P.child2 = child2;
            P.idx1 = idx1;
            P.idx2 = idx2;

        end

    end

end
