classdef surfacefun
%SURFACEFUN   A class for representing functions on surfaces.
%   SURFACEFUN(DOM) constructs a SURFACEFUN representing the zero function
%   over the SURFACEMESH DOM.
%
%   SURFACEFUN(FUNC, DOM) constructs a SURFACEFUN representing the function
%   handle FUNC over the SURFACEMESH DOM. The function handle must be of
%   the form @(x,y,z) ..., i.e., a function of three-dimensional Cartiesian
%   coordinates.
%
%   SURFACEFUN(VALS, DOM) constructs a SURFACEFUN from a cell array of
%   function values at tensor-product Chebyshev nodes.

    properties

        domain % Domain respresenting the surface patches
        vals   % Functions values

    end

    methods

        function obj = surfacefun(varargin)

            func = @(x,y,z) 0*x;
            vals = {};

            isValidDom = @(dom) isa(dom, 'surfacemesh');
            if ( nargin == 0 )
                % Empty SURFACEFUN:
                obj.vals = {};
                return
            elseif ( isa(varargin{1}, 'surfacefun') )
                error('Invalid call to surfacefun constructor.');
            elseif ( nargin == 1 && isValidDom(varargin{1}) )
                % Call is: SURFACEFUN(DOM)
                dom = varargin{1};
            elseif ( nargin == 2 && isValidDom(varargin{2}) )
                dom = varargin{2};
                if ( isa(varargin{1}, 'function_handle') )
                    % Call is: SURFACEFUN(FUNC, DOM)
                    func = varargin{1};
                elseif ( iscell(varargin{1}) && length(varargin{1}) == length(dom) )
                    % Call is: SURFACEFUN(VALS, DOM)
                    vals = varargin{1};
                    if ( size(vals, 2) > 1 )
                        % Return an array of surfacefuns:
                        obj = [];
                        for k = 1:size(vals, 2)
                            obj = [obj surfacefun(vals(:,k), dom)]; %#ok<AGROW>
                        end
                        return
                    end
                elseif ( isscalar(varargin{1}) )
                    vals = cell(length(dom), 1);
                    for k = 1:length(dom)
                        vals{k} = repmat(varargin{1}, size(dom.x{k}));
                    end
                else
                    error('SURFACEFUN:surfacefun:invalid', ...
                        'Invalid call to surfacefun constructor.');
                end
            else
                error('SURFACEFUN:surfacefun:invalid', ...
                    'Invalid call to SURFACEFUN constructor.');
            end

            obj.domain = dom;

            if ( isempty(vals) )
                n = order(dom)+1;
                xx = reshape([dom.x{:}], [], 1);
                yy = reshape([dom.y{:}], [], 1);
                zz = reshape([dom.z{:}], [], 1);
                ff = feval(func, xx, yy, zz);
                ff = reshape(ff, [n n length(dom)]);
                vals = cell(length(dom), 1);
                for k = 1:length(dom)
                    vals{k} = ff(:,:,k);
                end
            end
            obj.vals = vals;

        end

        function n = numArgumentsFromSubscript(obj,s,indexingContext) %#ok<INUSD>
        %NUMARGUMENTSFROMSUBSCRIPT   Number of arguments for customized indexing methods.
        %   Overloading NUMEL() gives the wrong NARGOUT for SUBSREF().
        %   Defining this function fixes it.
        %
        % See also NUMEL, NARGOUT, SUBSREF.
            n = 1;
        end

    end

    methods ( Static )

        C = vals2coeffs(V);
        V = coeffs2vals(C);

    end

end
