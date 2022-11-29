function K = krons(varargin)

if ( nargin == 1 )
    K = varargin{:};
else
    K = kron(varargin{1}, krons(varargin{2:end}));
end

end
