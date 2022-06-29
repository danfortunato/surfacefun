function normF = norm(f, varargin)
%NORM   Norm of a SURFACEFUN.
%   For SURFACEFUN objects:
%       NORM(F) = sqrt(integral of abs(F)^2).
%       NORM(F, 2) is the same as NORM(F).
%       NORM(F, 'fro') is also the same as NORM(F).
%       NORM(F, 1) = integral of abs(F).
%       NORM(F, P) = (integral of abs(F)^P)^(1/P).
%       NORM(F, inf) = estimated global maximum in absolute value.
%       NORM(F, -inf) = estimated global minimum in absolute value.
%       NORM(F, 'max') is the same as NORM(F, inf).
%       NORM(F, 'min') is the same as NORM(F, -inf).
%       NORM(F, 'H1') = sqrt( ||f||_2^2 + ||grad(f)||_2^2 )
%       NORM(F, 'H2') = sqrt( ||f||_H1^2 + ||grad(diffx(f))||_2^2
%                                        + ||grad(diffy(f))||_2^2
%                                        + ||grad(diffz(f))||_2^2 )
%       NORM(F, 'lap') = sqrt( ||f||_2 + ||lap(f)||_2^2 )
%
%   NORM(F, 'all') and NORM(F, P, 'all') return an array of norms on each
%   patch of F.

% Parse arguments.
p = 2;
reduce = true;
if ( nargin == 2 )
    if ( strcmp(varargin{1}, 'all') )
        reduce = false;
    else
        p = varargin{1};
    end
elseif ( nargin == 3 )
    p = varargin{1};
    if ( strcmp(varargin{2}, 'all') )
        reduce = false;
    end
end

% Empty SURFACEFUN has norm 0.
if ( isempty(f) )
    normF = 0;
    return
end

switch ( p )
    case 'fro'
        normF = norm(f, 2);
        reduceFun = @(x) sqrt( sum(x.^2) );

    case {inf, 'inf', 'max'}
        normF = cellfun(@(u) max(abs(u(:))), f.vals);
        reduceFun = @max;

    case {-inf, '-inf', 'min'}
        normF = cellfun(@(u) min(abs(u(:))), f.vals);
        reduceFun = @min;

    case 'H1'
        [fx, fy, fz] = grad(f);
        norm_L2 = norm(f,  'all');
        norm_x  = norm(fx, 'all');
        norm_y  = norm(fy, 'all');
        norm_z  = norm(fz, 'all');
        normF = sqrt( norm_L2.^2 + norm_x.^2 + norm_y.^2 + norm_z.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    case 'H2'
        [fx,  fy,  fz]  = grad(f);
        [fxx, fxy, fxz] = grad(fx);
        [fyx, fyy, fyz] = grad(fy);
        [fzx, fzy, fzz] = grad(fz);
        norm_H1 = norm(f, 'H1', 'all');
        norm_xx = norm(fxx, 'all'); norm_xy  = norm(fxy, 'all'); norm_xz  = norm(fxz, 'all');
        norm_yx = norm(fyx, 'all'); norm_yy  = norm(fyy, 'all'); norm_yz  = norm(fyz, 'all');
        norm_zx = norm(fzx, 'all'); norm_zy  = norm(fzy, 'all'); norm_zz  = norm(fzz, 'all');
        normF = sqrt( norm_H1.^2 + norm_xx.^2 + norm_xy.^2 + norm_xz.^2 + ...
                                   norm_yx.^2 + norm_yy.^2 + norm_yz.^2 + ...
                                   norm_zx.^2 + norm_zy.^2 + norm_zz.^2);
        reduceFun = @(x) sqrt( sum(x.^2) );

    case 'lap'
        norm_L2  = norm(f,      'all');
        norm_lap = norm(lap(f), 'all');
        normF = sqrt( norm_L2.^2 + norm_lap.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    otherwise
        if ( isnumeric(p) && isreal(p) )
            if ( abs(round(p) - p) < eps )
                normF = integrate(f, p);
                reduceFun = @(x) ( sum(x.^p) ).^(1/p);
            else
                error('SURFACEFUN:norm:norm', ...
                    'SURFACEFUN does not support this norm.');
            end
        else
            error('SURFACEFUN:norm:unknown', 'Unknown norm.');
        end

end

% Combine norms on each element.
if ( reduce )
    normF = reduceFun(normF);
end

end

function int = integrate(f, p)
%INTEGRATE   Compute (integral of abs(F)^P)^(1/P)

% P should be an integer.
p = round(p);

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N*P points.
int = zeros(length(f), 1);
for k = 1:length(f)
    [nv, nu] = size(f.vals{k});
    qu = nu*p;
    qv = nv*p;
    nv = min(nv, qv);
    nu = min(nu, qu);
    coeffs = chebtech2.vals2coeffs(chebtech2.vals2coeffs(f.vals{k}).').';
    U = zeros(qv, qu);
    U(1:nv,1:nu) = coeffs(1:nv,1:nu);
    V = chebtech2.coeffs2vals(chebtech2.coeffs2vals(U).').';
    J_cfs = chebtech2.vals2coeffs(chebtech2.vals2coeffs(f.domain.J{k}).').';
    J_q = zeros(qv, qu);
    J_q(1:nv,1:nu) = J_cfs(1:nv,1:nu);
    J = chebtech2.coeffs2vals(chebtech2.coeffs2vals(J_q).').';
    wu = chebtech2.quadwts(qu); wu = wu(:);
    wv = chebtech2.quadwts(qv); wv = wv(:);
    int(k) = sum(sum(abs(V).^p .* wv .* wu.' .* sqrt(J))).^(1/p);
end

end
