function f = uminus(f)
%-   Unary minus for a SURFACEFUNV.
%   -F negates the SURFACEFUNV F.
%
%   G = UMINUS(F) is called for the syntax '-F'.
%
%   See also UPLUS.

if ( isempty(f) )
    return
end

f.components = cellfun(@uminus, f.components, 'UniformOutput', false);

end
