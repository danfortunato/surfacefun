classdef surfacemesh < handle
%SURFACEMESH

    properties

        x
        y
        z

    end

    properties ( Hidden )

        xu
        yu
        zu
        xv
        yv
        zv

        ux
        uy
        uz
        vx
        vy
        vz

        E
        F
        G
        J
        singular

        facenormals
        edgenormals
        connectivity

    end

    methods

        function dom = surfacemesh(x, y, z)

            if ( nargin == 0 )
                return
            end

            if ( ~(iscell(x) && iscell(y) && iscell(z) && ...
                   all(size(x) == size(y)) && all(size(y) == size(z))) )
                error('X, Y, and Z must be cell arrays of the same size.');
            end

            dom.x = x;
            dom.y = y;
            dom.z = z;

            xu = cell(size(x)); xv = cell(size(x));
            yu = cell(size(x)); yv = cell(size(x));
            zu = cell(size(x)); zv = cell(size(x));
            ux = cell(size(x)); vx = cell(size(x));
            uy = cell(size(x)); vy = cell(size(x));
            uz = cell(size(x)); vz = cell(size(x));
            E = cell(size(x));
            F = cell(size(x));
            G = cell(size(x));
            J = cell(size(x));

            n = size(x{1}, 1);
            D = diffmat(n);
            singular = false(size(x));
            for k = 1:length(x)
                xu{k} = x{k} * D.'; xv{k} = D * x{k};
                yu{k} = y{k} * D.'; yv{k} = D * y{k};
                zu{k} = z{k} * D.'; zv{k} = D * z{k};
                E{k} = xu{k}.*xu{k} + yu{k}.*yu{k} + zu{k}.*zu{k};
                G{k} = xv{k}.*xv{k} + yv{k}.*yv{k} + zv{k}.*zv{k};
                F{k} = xu{k}.*xv{k} + yu{k}.*yv{k} + zu{k}.*zv{k};
                J{k} = E{k}.*G{k} - F{k}.^2;

                scl = max(abs(G{k}.*xu{k}-F{k}.*xv{k}), [], 'all');
                if ( any(abs(J{k}) < 1e-10*scl, 'all') )
                    singular(k) = true;
                    ux{k} = G{k}.*xu{k}-F{k}.*xv{k};
                    uy{k} = G{k}.*yu{k}-F{k}.*yv{k};
                    uz{k} = G{k}.*zu{k}-F{k}.*zv{k};
                    vx{k} = E{k}.*xv{k}-F{k}.*xu{k};
                    vy{k} = E{k}.*yv{k}-F{k}.*yu{k};
                    vz{k} = E{k}.*zv{k}-F{k}.*zu{k};
                else
                    ux{k} = (G{k}.*xu{k}-F{k}.*xv{k})./J{k};
                    uy{k} = (G{k}.*yu{k}-F{k}.*yv{k})./J{k};
                    uz{k} = (G{k}.*zu{k}-F{k}.*zv{k})./J{k};
                    vx{k} = (E{k}.*xv{k}-F{k}.*xu{k})./J{k};
                    vy{k} = (E{k}.*yv{k}-F{k}.*yu{k})./J{k};
                    vz{k} = (E{k}.*zv{k}-F{k}.*zu{k})./J{k};
                end
            end

            dom.xu = xu; dom.xv = xv; dom.ux = ux; dom.vx = vx;
            dom.yu = yu; dom.yv = yv; dom.uy = uy; dom.vy = vy;
            dom.zu = zu; dom.zv = zv; dom.uz = uz; dom.vz = vz;
            dom.E = E;
            dom.F = F;
            dom.G = G;
            dom.J = J;
            dom.singular = singular;

            dom.connectivity = buildConnectivity(dom);
            dom.facenormals = computeNormals(dom);
            %dom.edgenormals = normal(dom, 'edges');

        end

        function n = numArgumentsFromSubscript(varargin)
        %NUMARGUMENTSFROMSUBSCRIPT   Number of arguments for customized indexing methods.
        %   Overloading NUMEL() gives the wrong NARGOUT for SUBSREF().
        %   Defining this function fixes it.
            n = 1;
        end

    end

    methods ( Static )

        dom = sphere(varargin);
        dom = ellipsoid(varargin);
        dom = hemisphere(varargin);
        dom = blob(varargin);
        dom = torus(varargin);
        dom = sharptorus(varargin);
        dom = cyclide(varargin);
        dom = stellarator(varargin);
        dom = cube(varargin);
        dom = flat_sphere(varargin);
        dom = teardrop(varargin);
        dom = mobius(varargin);
        dom = fromGmsh(varargin);
        dom = fromRhino(varargin);

    end

end

% function [Dx, Dy, Dz] = diffs(x, y, z)
% %DIFFS   Differential operators on surfaces.
% %   [DX, DY, DZ] = DIFFS(X, Y, Z) returns the tangential derivative
% %   operators DX, DY, and DZ over the surface defined by the coordinates
% %   (X, Y, Z). The coordinates must be given as N x N matrices sampled at
% %   tensor-product Chebyshev nodes. DX, DY, and DZ are operators of size
% %   N^2 x N^2.
% 
% persistent D
% n = size(x, 1);
% if ( size(D, 1) ~= n )
%     D = diffmat(n);
% end
% 
% xu = x * D.'; xv = D * x;
% yu = y * D.'; yv = D * y;
% zu = z * D.'; zv = D * z;
% E = xu.*xu + yu.*yu + zu.*zu;
% G = xv.*xv + yv.*yv + zv.*zv;
% F = xu.*xv + yu.*yv + zu.*zv;
% J = E.*G - F.^2;
% ux = (G.*xu-F.*xv)./J; vx = (E.*xv-F.*xu)./J;
% uy = (G.*yu-F.*yv)./J; vy = (E.*yv-F.*yu)./J;
% uz = (G.*zu-F.*zv)./J; vz = (E.*zv-F.*zu)./J;
% I = eye(n); % or speye?
% Du = kron(D, I);
% Dv = kron(I, D);
% Dx = ux(:).*Du + vx(:).*Dv;
% Dy = uy(:).*Du + vy(:).*Dv;
% Dz = uz(:).*Du + vz(:).*Dv;
% 
% end
