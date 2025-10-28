function dom = import(filename, format, varargin)
%IMPORT   Import a surfacemesh.

if ( ~exist(filename, 'file') )
    error('File not found: ''%s''', filename);
end

if ( nargin == 1 && (ischar(filename) || isstring(filename)) )
    [~, ~, ext] = fileparts(filename);
    switch ext
        case '.msh'
            format = 'gmsh';
        case '.csv'
            format = 'rhino';
        otherwise
            error('Cannot determine mesh format from file extension.');
    end
end

switch lower(format)
    case 'gmsh'
        dom = surfacemesh.import.gmsh(filename, varargin{:});
    case {'csv', 'rhino', 'rhinoceros'}
        dom = surfacemesh.import.csv(filename, varargin{:});
    otherwise
        error('Unknown mesh format ''%s''.', format);
end

end
