function f = ctranspose(f)
%'   Complex conjugate transpose of a SURFACEFUN.
%   F' is the complex conjugate transpose of the SURFACEFUN F.

% Empty check:
if ( isempty(f) )
    return
end

f = conj(f).';

end
