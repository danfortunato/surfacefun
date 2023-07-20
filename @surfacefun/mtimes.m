function f = mtimes(f, g)
%*   Matrix multiplication for SURFACEFUN.
%   c*F or F*c multiplies the SURFACEFUN F by the scalar c. F may be an
%   array-valued SURFACEFUN and c may be a matrix of scalars with
%   compatible dimensions.
%
%   See also TIMES.

if ( isempty(f) || isempty(g) )
    f = [];
elseif ( isa(f, 'surfacefun') && isnumeric(g) )
    if ( isscalar(g) )
        f = times(f, g);
    else
        % Make sure f is a quasimatrix:
        if ( size(f, 1) > 1 && size(f, 2) > 1 )
            error('SURFACEFUN:mtimes:notquasimatrix', 'Not a quasimatrix.');
        end

        % Check the dimensions:
        [m, n] = size(g);
        if ( m ~= size(f, 2) || ~ismatrix(g) )
            error('SURFACEFUN:mtimes:dims', 'Matrix dimensions must agree. Did you mean ''.*''?');
        end

        out = repmat(0*f(:,1), 1, size(g, 2));
        for k = 1:n
            for j = 1:m
                out(:,k) = out(:,k) + f(:,j).*g(j,k);
            end
        end
        f = out;
    end
elseif ( isnumeric(f) && isa(g, 'surfacefun') )
    f = mtimes(g.', f.').';
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefun') )
    % MTIMES will take the L^2 inner product between arrays of surfacefuns
    if ( ismatrix(f) && ismatrix(g) )
        % Check the dimensions:
        if ( size(f, 2) ~= size(g, 1) )
            error('SURFACEFUN:mtimes:dims', 'Matrix dimensions must agree. Did you mean ''.*''?');
        end
        out = zeros(size(f, 1), size(g, 2));
        for i = 1:size(f, 1)
            for j = 1:size(g, 2)
                for k = 1:size(f, 2)
                    out(i,j) = out(i,j) + integral2( f(i,k).*conj(g(k,j)) );
                end
            end
        end
        f = out;
    end
else
    error('SURFACEFUN:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end

end
