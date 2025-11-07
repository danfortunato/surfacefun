function initialize(S, varargin)
%INITIALIZE   Initialize and solve local subproblems.
%   INITIALIZE(S, RHS) will initialize the SURFACEOP object S on each of
%   the subpatches of S.domain using the righthand side RHS.
%
%   INITIALIZE(S) assumes the problem is homogeneous (i.e., RHS = 0).
%
%   The full sequence for solving a problem using a SURFACEOP object S
%   is:
%
%      initialize(S, RHS)
%      build(S)   % (optional)
%      sol = S\bc % or sol = solve(S, bc)
%
%   See also BUILD, SOLVE.

% Initialize all leaf patches:
switch S.method
    case 'DtN'
        S.patches = surfaceop.leaf.initialize_DtN(S.op, S.domain, varargin{:});
    case 'ItI'
        S.patches = surfaceop.leaf.initialize_ItI(S.op, S.domain, S.eta, varargin{:});
end

% Handle hanging nodes by interpolating any coarse leaf boundary operators to
% the split edges specified in S.domain.split. S.domain.split{k}(j) is true if
% the j-th side of the k-th element is split, with sides stored in the order:
% left, right, down, up. For example, for the mesh:
%
%    +-------+---+
%    |       | 2 |
%    |   1   +---+
%    |       | 3 |
%    +-------+---+
%
% we would have:
%
%    S.domain.split{1} = [0 1 0 0];
%    S.domain.split{2} = [0 0 0 0];
%    S.domain.split{3} = [0 0 0 0];
%
% Note: This will break if the mesh is not 2:1 level restricted.
S.patches = splitHangingEdges(S.patches, S.domain);

end

function patches = splitHangingEdges(patches, dom)

if ( isempty(patches) )
    return
end

n = patches{1}.n;
nskel = n-2;
I_skel = eye(nskel);
[xskel,      ~, vskel]      = chebpts(nskel, [-1 1], 1);
[xskel_lsub, ~, vskel_lsub] = chebpts(nskel, [-1 0], 1);
[xskel_rsub, ~, vskel_rsub] = chebpts(nskel, [ 0 1], 1);
w = chebtech1.quadwts(nskel);
w = w(:);

% Interpolate to two type-1 grids on [-1 0 1] from one type-1 grid on [-1 1]
B_2f1 = barymat([xskel_lsub ; xskel_rsub], xskel, vskel);

% Interpolate to one type-1 grid on [-1 1] from two type-1 grids on [-1 0 1]
% Note: Maybe we should average the endpoints here?
B_1f2l = barymat(xskel(xskel<=0), xskel_lsub, vskel_lsub);
B_1f2r = barymat(xskel(xskel>0),  xskel_rsub, vskel_rsub);
B_1f2  = blkdiag(B_1f2l, B_1f2r);

% Interpolate to the midpoint from one type-2 grid on [-1 1]
[xleaf, ~, wleaf] = chebpts(n, 2);
B_mid = barymat(0, xleaf, wleaf);

for k = 1:length(patches)
    split = dom.split{k};
    if ( any(split) )
        % We have a hanging node.
        % Map each split side to two grids:
        B_split_1f2 = repmat({I_skel}, 4, 1);
        B_split_2f1 = repmat({I_skel}, 4, 1);
        B_split_1f2(split) = {B_1f2};
        B_split_2f1(split) = {B_2f1};
        B_split_1f2 = blkdiag(B_split_1f2{:});
        B_split_2f1 = blkdiag(B_split_2f1{:});
        patches{k}.S        = patches{k}.S * B_split_1f2;
        patches{k}.BtB      = B_split_2f1 * patches{k}.BtB * B_split_1f2;
        patches{k}.du_part  = B_split_2f1 * patches{k}.du_part;
        patches{k}.normal_d = B_split_2f1 * patches{k}.normal_d;
        patches{k}.xyz      = B_split_2f1 * patches{k}.xyz;
        nsides = sum(split+1);
        patches{k}.BtB_scl = repmat({ones(nskel,1)}, nsides, 1);

        % Insert edge vertices for the split edges:
        x = dom.x{k};
        y = dom.y{k};
        z = dom.z{k};
        edges = patches{k}.edges;
        edges_left = edges(1,:);
        w_left = w;
        if ( split(1) )
            mid = B_mid * [x(:,1) y(:,1) z(:,1)];
            edges_left = [ x(1,1) y(1,1) z(1,1) mid(1) mid(2) mid(3) nskel ;
                           mid(1) mid(2) mid(3) x(n,1) y(n,1) z(n,1) nskel ];
            w_left = [ 0.5*w ; 0.5*w ];
        end
        edges_right = edges(2,:);
        w_right = w;
        if ( split(2) )
            mid = B_mid * [x(:,n) y(:,n) z(:,n)];
            edges_right = [ x(1,n) y(1,n) z(1,n) mid(1) mid(2) mid(3) nskel ;
                            mid(1) mid(2) mid(3) x(n,n) y(n,n) z(n,n) nskel ];
            w_right = [ 0.5*w ; 0.5*w ];
        end
        edges_down = edges(3,:);
        w_down = w;
        if ( split(3) )
            mid = B_mid * [x(1,:).' y(1,:).' z(1,:).'];
            edges_down = [ x(1,1) y(1,1) z(1,1) mid(1) mid(2) mid(3) nskel ;
                           mid(1) mid(2) mid(3) x(1,n) y(1,n) z(1,n) nskel ];
            w_down = [ 0.5*w ; 0.5*w ];
        end
        edges_up = edges(4,:);
        w_up = w;
        if ( split(4) )
            mid = B_mid * [x(n,:).' y(n,:).' z(n,:).'];
            edges_up = [ x(n,1) y(n,1) z(n,1) mid(1) mid(2) mid(3) nskel ;
                         mid(1) mid(2) mid(3) x(n,n) y(n,n) z(n,n) nskel ];
            w_up = [ 0.5*w ; 0.5*w ];
        end
        patches{k}.edges = [ edges_left ; edges_right ; edges_down ; edges_up ];
        wskel = [ w_left ; w_right ; w_down ; w_up ];
        JJ = patches{k}.w ./ [w ; w ; w ; w];
        patches{k}.w = wskel .* (B_split_2f1 * JJ);
    end
end

end
