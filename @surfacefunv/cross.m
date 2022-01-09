function h = cross(f, g, varargin)
%CROSS   Vector cross product.
%   CROSS(F, G) returns a SURFACEFUNV representing the 3D cross product of
%   the vector fields F and G, where F and G are SURFACEFUNV objects or
%   constant vectors.
%
%   See also DOT.

% Support for CROSS(F, G, H, ...)
if ( nargin > 2 )
    h = cross(f, cross(g, varargin{:}));
    return
end

% Empty check:
if ( isempty(f) || isempty(g) )
    h = surfacefunv;
    return
end

if ( isa(f, 'surfacefunv') && isa(g, 'surfacefunv') )
    % Cross product of two SURFACEFUNVs:
    fc = f.components;
    gc = g.components;
    h = surfacefunv(fc{2}.*gc{3} - fc{3}.*gc{2}, ...
                    fc{3}.*gc{1} - fc{1}.*gc{3}, ...
                    fc{1}.*gc{2} - fc{2}.*gc{1});
elseif ( isa(f, 'surfacefunv') && isnumeric(g) && numel(g) == 3 )
    % Cross SURFACEFUNV with vector:
    fc = f.components;
    h = surfacefunv(fc{2}.*g(3) - fc{3}.*g(2), ...
                    fc{3}.*g(1) - fc{1}.*g(3), ...
                    fc{1}.*g(2) - fc{2}.*g(1));
elseif ( isnumeric(f) && numel(f) == 3  && isa(g, 'surfacefunv') )
    % Cross vector with SURFACEFUNV:
    gc = g.components;
    h = surfacefunv(f(2).*gc{3} - f(3).*gc{2}, ...
                    f(3).*gc{1} - f(1).*gc{3}, ...
                    f(1).*gc{2} - f(2).*gc{1});
else
    error('SURFACEFUNV:cross:invalid', ...
        'F and G must be surfacefunv objects or constant vectors.')
end

end
