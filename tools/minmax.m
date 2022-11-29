function varargout = minmax(X)
%MINMAX   Minimum and maximum elements of an array.
%   MINMAX(X) returns a 1x2 row vector containing the minimum and maximum
%   elements of X.
%
%   [XMIN, XMAX] = MINMAX(X) instead returns the minimum and maximum
%   elements of X separately.


xmin = min(X, [], 'all');
xmax = max(X, [], 'all');

if ( nargout < 2 )
    varargout{1} = [xmin xmax];
else
    varargout{1} = xmin;
    varargout{2} = xmax;
end

end
