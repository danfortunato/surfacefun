function varargout = subsref(f, index)
%SUBSREF   Subscripted reference for a SURFACEFUN.
%   F(X, Y, Z) returns the values of the SURFACEFUN F evaluated at (X, Y,
%   Z). See FEVAL for further details.
%
%   F.PROP returns the property PROP of the SURFACEFUN F as defined by
%   GET(SOL, 'PROP').
%
%   See also FEVAL, GET.

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'

        if ( all(size(f) ~= 1) )
            % Where to evaluate:
            x = idx{1};
            y = idx{2};
            z = idx{3};
            out = feval(f, x, y, z);
        else
            [varargout{1:nargout}] = builtin('subsref', f, index);
            return
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '.'

        % Call GET() for .PROP access.
        out = get(f, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end

    otherwise

        error('SURFACEFUN:subsref:unexpectedType',...
            ['Unexpected index.type of ', index(1).type]);
end

% Convert to a cell:
varargout = {out};

end
