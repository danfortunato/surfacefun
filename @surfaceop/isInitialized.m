function out = isInitialized(S)
%ISINITIALIZED   Check to see if a SURFACEOP has been initialized.

out = ~isempty(S.patches) && (size(S.patches{1}.S, 2) > 0 || ~isempty(S.patches{1}.u_part) );

end
