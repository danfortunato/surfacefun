function f = subsasgn(f, index, val)
%SUBSASGN   Subscripted assignment to a SURFACEFUN object.
%   F.PROP = VAL allows property assignment for F.DOMAIN and F.VALS.
%
%   F.VEC(:) = X sets the function values of the SURFACEFUN F from the
%   column vector X of length NUMEL(F). If F is an array-valued SURFACEFUN
%   with M columns, then X should be a matrix of size NUMEL(F(1)) x M.
%
%   See also VEC.

switch index(1).type
    case '.'
        if ( strcmpi(index(1).subs, 'vec') )
            [nu, nv] = size(f(1).domain);
            npts = nu*nv;
            numPatches = length(f(1).domain);
            numDOF = numel(f(1));
            numFuns = size(f, 2);
            if ( ~all(size(val) == [numDOF numFuns]) )
                error('SURFACEFUN:subsasgn:dimensions', ...
                    'Data and surfacefun dimensions are incompatible.');
            end
            for j = 1:numFuns
                for k = 1:numPatches
                    f(:,j).vals{k}(:) = val((k-1)*npts+(1:npts), j);
                end
            end
        else
            f = builtin('subsasgn', f, index, val);
        end

    case '()'
        f = builtin('subsasgn', f, index, val);

    otherwise
        error('SURFACEFUN:subsasgn:unexpectedType', ...
            ['Unexpected index.type of ' index(1).type]);
end

end
