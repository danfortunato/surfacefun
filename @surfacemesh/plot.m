function varargout = plot(dom, varargin)
%PLOT   Plot a surface.
%   PLOT(DOM) gives a surface plot of the SURFACEMESH DOM, the same as
%   WIREFRAME(DOM).
%
%   See also WIREFRAME, MESH.

[varargout{1:nargout}] = wireframe(dom, varargin{:});

end
