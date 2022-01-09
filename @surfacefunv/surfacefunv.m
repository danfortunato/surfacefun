classdef surfacefunv
%SURFACEFUNV   A class for representing functions on surfaces.
%   SURFACEFUNV(DOM) constructs an SURFACEFUNV over the SURFACEMESH DOM.
%
%   SURFACEFUNV(F, G, H) constructs a SURFACEFUNV representing the vector
%   field (F, G, H), where F, G, and H are SURFACEFUNs with the same
%   domain.
%
%   SURFACEFUNV(F, G, H, DOM) constructs a SURFACEFUNV representing the
%   vector field (F, G, H) over the SURFACEMESH DOM, where F, G, and H are
%   function handles of the form @(x,y,z) ... or cell arrays of function
%   values at tensor-product Chebyshev nodes.

    properties

        components
        isTransposed

    end

    methods

        function f = surfacefunv(varargin)

            if ( nargin == 0 )
                % Empty SURFACEFUNV:
                f.components = {surfacefun, surfacefun, surfacefun};
            elseif ( nargin == 1 && isa(varargin{1}, 'surfacemesh') )
                % Call is: SURFACEFUNV(DOM)
                dom = varargin{1};
                f.components = {surfacefun(dom), ...
                                surfacefun(dom), ...
                                surfacefun(dom)};
            elseif ( nargin == 3 && isa(varargin{1}, 'surfacefun') && ...
                                    isa(varargin{2}, 'surfacefun') && ...
                                    isa(varargin{3}, 'surfacefun'))
                % Call is: SURFACEFUNV(F, G, H)
                f.components = varargin;
            elseif ( nargin == 4 && isa(varargin{4}, 'surfacemesh') )
                % Call is: SURFACEFUNV(F, G, H, DOM)
                dom = varargin{4};
                f.components = {surfacefun(varargin{1}, dom), ...
                                surfacefun(varargin{2}, dom), ...
                                surfacefun(varargin{3}, dom)};
            else
                error('SURFACEFUNV:surfacefunv:invalid', ...
                        'Invalid call to surfacefunv constructor.');
            end

            f.isTransposed = false;

        end

    end

end
