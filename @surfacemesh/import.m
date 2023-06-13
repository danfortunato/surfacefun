function dom = import(filename, format, varargin)
%IMPORT   Import a surfacemesh.

switch lower(format)
    case 'gmsh'
        dom = surfacemesh.fromGmsh(filename, varargin{:});
    case {'rhino', 'rhinoceros'}
        dom = surfacemesh.fromRhino(filename, varargin{:});
    otherwise
        error('Unknown mesh format ''%s''.', format);
end

end
