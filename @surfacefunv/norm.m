function normF = norm(f)
%NORM   Frobenius norm of a SURFACEFUNV.
%   NORM(F) = sqrt(Fx.^2 + Fy.^2 + Fz.^2), where F is a SURFACEFUNV with
%   components Fx, Fy, and Fz.

% Empty check:
if ( isempty(f) )
    normF = [];
    return
end

fc = f.components;
normF = sqrt(fc{1}.^2 + fc{2}.^2 + fc{3}.^2);

end
