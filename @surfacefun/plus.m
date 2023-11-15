function h = plus(f, g)
%+   Plus for SURFACEFUN.
%   F + G adds the SURFACEFUN F and G. F and G must have the same domains
%   and discretization sizes. F and G may also be scalars.
%
%   See also MINUS.

% Empty check:
if ( isempty(f) )
    return
end

h = compose(@plus, f, g);

end
