function varargout = mldivide(varargin)
%MLDIVIDE   Perform a strip solve.
%   MLDIVIDE is a shorthand for SOLVE().
%
%   See also SOLVE.

[varargout{1:nargout}] = solve(varargin{:});

end
