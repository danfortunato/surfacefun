function chain(varargin)

for k = 1:nargin
    varargin{k}();
end

end
