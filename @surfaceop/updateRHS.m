function S = updateRHS(S, f)
%UPDATERHS   Update the RHS function in a SURFACEOP.
%   UPDATERHS(S, F) or S.RHS = F updates the RHS function F in the
%   SURFACEOP S. The function F may be a constant, a function handle, a
%   SURFACEFUN object defined on S, or a cell array containing tensor-
%   product Chebyshev values for each patch. This is useful as only a small
%   fraction of the initialization phase needs to be repeated.
%   See also BUILD.

assert(isInitialized(S), 'The surfaceop has not been initialized.')

% Extract values if f is a surfacefun:
if ( isa(f, 'surfacefun') )
    if ( size(f, 2) == 1 )
        f = f.vals;
    else
        vals = cell(length(f.domain), size(f, 2));
        for k = 1:size(f, 2)
            vals(:,k) = f(:,k).vals;
        end
        f = vals;
    end
end

% Duplicate f if it is a scalar or function handle:
if ( ~iscell(f) )
    f = repmat({f}, length(S.domain), 1);
end

switch S.method
    case 'DtN'
        updateRHS = @updateRHS_DtN;
    case 'ItI'
        updateRHS = @updateRHS_ItI;
end

% Update RHS of each patch:
if ( isBuilt(S) )
    S.patches{1} = updateRHS(S.patches{1}, f);
else
    for k = 1:numel(S.patches)
        S.patches{k} = updateRHS(S.patches{k}, f(k,:));
    end
end

end
