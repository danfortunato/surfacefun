function f = uminus(f)
%-   Unary minus for a SURFACEFUN.
%   -F negates the SURFACEFUN F.
%
%   G = uminus(F) is called for the syntax '-F'.
%
%   See also UPLUS.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@uminus, f);

end
