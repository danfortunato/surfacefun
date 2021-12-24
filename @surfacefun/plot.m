function varargout = plot(varargin)
%PLOT   Plot a SURFACEFUN.
%   PLOT(F) gives a surface plot of the SURFACEFUN F, the same as SURF(F).
%
%   See also SURF, MESH, CONTOUR.

[varargout{1:nargout}] = surf(varargin{:});

end
