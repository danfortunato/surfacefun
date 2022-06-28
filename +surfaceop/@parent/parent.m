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
        flip1
        flip2
        scl1
        scl2
        A           % Interface linear system
        dA          % Decomposition of interface linear system

    end

    methods

        function P = parent(domain, id, S, D2N, D2N_scl, u_part, du_part, A, dA, edges, xyz, child1, child2, idx1, idx2, flip1, flip2, scl1, scl2)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = domain;
            P.id = id;
            P.S = S;
            P.D2N = D2N;
            P.D2N_scl = D2N_scl;
            P.u_part = u_part;
            P.du_part = du_part;
            P.A = A;
            P.dA = dA;
            P.edges = edges;
            P.xyz = xyz;
            P.len = child1.len + child2.len;

            % Assign children:
            P.child1 = child1;
            P.child2 = child2;
            P.idx1 = idx1;
            P.idx2 = idx2;
            P.flip1 = flip1;
            P.flip2 = flip2;
            P.scl1 = scl1;
            P.scl2 = scl2;

        end

    end

end
