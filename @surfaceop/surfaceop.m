classdef surfaceop < handle
%SURFACEOP   An SEM solver on a surface.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain = []   % Surface mesh.
        op            % Differential operator.
        patches = {}  % Cell array of surfaceop.Patches.
        mergeIdx = {} % Cell array of merge indices.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = surfaceop(dom, op, rhs)
        %SURACEFOP   Constructor for the SURFACEOP class.

            if ( nargin == 0 )
                % Construct empty object:
                return
            end

            if ( nargin < 3 )
                rhs = 0;
            end

            obj.op = parsePDO(op);
 
            % Assign the domain:
            obj.domain = dom;

            % Build merge indices:
            obj.mergeIdx = surfaceop.defaultIdx(dom);

            % Initialize patches:
            obj.initialize(rhs);

        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Static )
        
        function mergeIdx = defaultIdx(dom)
        %DEFAULTIDX   Make a default (naive) merge index for a domain.
        %   MERGEIDX = DEFAULTIDX(DOM) automatically constructs a merge index
        %   MERGEIDX for the domain DOM. The merge index will take the form
        %
        %       MERGEIDX = {[1 2 ; 3 4 ; ... ; N-1 N],
        %                   [1 2 ; ... ; N/2-1 N/2],
        %                    ...,
        %                   [1 2]},
        %
        %   where N is the number of patches in DOM. NaNs are introduced where
        %   there are an odd number of patches on a level.
        %
        %   MERGEIDX = DEFAULTIDX(N) is equivalent.

            % Parse input:
            if ( isnumeric(dom) && isscalar(dom) )
                np = dom;            % Number of patches.
            else
                np = length(dom);    % Number of patches.
            end
            nl = ceil(log2(np));      % Number of layers.
            mergeIdx = cell(1, nl);   % Initialize.

            for l = 1:nl
                if ( mod( np, 2 ) )
                    ii = [1:np, NaN]; % Pad with a NaN.
                else
                    ii = 1:np;
                end
                np = ceil( np/2 );    % Number of patches is halved each layer.
                mergeIdx{1,l} = reshape(ii, 2, np).';
            end

        end

    end

end

function out = parsePDO(in)

if ( ~isstruct(in) )
    error('Differential operator must be given as a struct.');
end

out = struct();
out.dxx = 0; out.dyy = 0; out.dzz = 0;
out.dxy = 0; out.dyz = 0; out.dxz = 0;
out.dx  = 0; out.dy  = 0; out.dz  = 0;
out.b = 0;

% Laplacian shorthand
if ( isfield(in, 'lap') )
    out.dxx = in.lap;
    out.dyy = in.lap;
    out.dzz = in.lap;
end

% Gradient shorthand
if ( isfield(in, 'grad') )
    out.dx = in.grad;
    out.dy = in.grad;
    out.dz = in.grad;
end

% Second derivatives
if ( isfield(in, 'dxx') ), out.dxx = in.dxx; end
if ( isfield(in, 'dyy') ), out.dyy = in.dyy; end
if ( isfield(in, 'dzz') ), out.dzz = in.dzz; end

% Cross derivatives
if ( isfield(in, 'dxy') ), out.dxy = in.dxy; end
if ( isfield(in, 'dyx') ), out.dxy = in.dyx; end
if ( isfield(in, 'dyz') ), out.dyz = in.dyz; end
if ( isfield(in, 'dzy') ), out.dyz = in.dzy; end
if ( isfield(in, 'dxz') ), out.dxz = in.dxz; end
if ( isfield(in, 'dzx') ), out.dxz = in.dzx; end

% First derivatives
if ( isfield(in, 'dx') ), out.dx = in.dx; end
if ( isfield(in, 'dy') ), out.dy = in.dy; end
if ( isfield(in, 'dz') ), out.dz = in.dz; end

% Zero-th derivative
if ( isfield(in, 'b') ), out.b = in.b; end

end
