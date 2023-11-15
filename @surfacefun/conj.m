function f = conj(f)
%CONJ   Comnplex conjugate of a SURFACEFUN.
%   CONJ(F) returns the complex conjugate of the SURFACEFUN F.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@conj, f);

end
