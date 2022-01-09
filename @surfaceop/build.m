function build(S)
%BUILD   Build the merge tree of patches.
%   BUILD(S) will construct the global solution operator from the local
%   operators in S.patches in the order defined in the merge tree
%   S.domain.mergeIdx. After building, S will have a single patch.
%
%   If S has not yet been initialized, then an error is thrown.
%
%   See also INITIALIZE, SOLVE.

if ( isa(S.patches, 'surfaceop') )
    % Recurse down and build lower level patches
    for k = 1:numel(S.patches)
        build(S.patches(k));
    end
    % Concatenate (since S already contains domain info.)
    S_sub = S.patches;
    S.patches = vertcat(S_sub.patches);
end

if ( ~isInitialized(S) )
    error('%s has not yet been initialized.', inputname(1));
end

% Build the patches:
%S.patches = build(S.domain, S.patches);

if ( length(S.domain) == 1 || numel(S.patches) == 1 )
    % Nothing to build for a single patch.
    return
end

% Loop over each level:
for j = 1:numel(S.mergeIdx) % <-- number of levels

    idxj = S.mergeIdx{j};               % Indicies for jth level.
    q = cell(size(idxj, 1),1);          % Initialize patches at jth level.

    % Perform all the merges at this level:
    for k = 1:size(idxj, 1)
        idxjk = idxj(k,:);              % Index for kth merge at jth level.
        idxjk(isnan(idxjk)) = [];       % Remove NaNs.
        pk = S.patches(idxjk);          % Extract patches.
        pk(cellfun(@isempty, pk)) = []; % Remove empty patches.
        if ( isempty(pk) )
            q{k} = [];                  % Empty merge.
        elseif ( j == numel(S.mergeIdx) )
            q{k} = merge(pk{:}, true);  % Perform top-level merge.
        else
            q{k} = merge(pk{:});        % Perform merge.
        end
    end

    S.patches = q; % Store back in S.patches for recursion.

end

end
