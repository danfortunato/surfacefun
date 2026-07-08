function idx = tuples(dim, n, type)
%TUPLES   Generate all tuples with a given sum.
%   IDX = TUPLES(DIM, N, '=') constructs a matrix IDX of all unique DIM-tuples
%   whose entries sum to N, in lexicographic order. IDX(:,j) is a vector of
%   length DIM. TUPLES(DIM, N) is equivalent to TUPLES(DIM, N, '=').
%
%   IDX = TUPLES(DIM, N, '<=') constructs a matrix IDX of all unique DIM-tuples
%   whose entries sum to at most N, in lexicographic order. IDX(:,j) is a vector
%   of length DIM.
%
%   TUPLES(DIM, N, '<') is equivalent to TUPLES(DIM, N-1, '<=').
%
%   Examples:
%
%   >> tuples(3, 2, '=')
%
%   ans =
% 
%        0     0     0     1     1     2
%        0     1     2     0     1     0
%        2     1     0     1     0     0
%
%   >> tuples(3, 2, '<=')
%
%   ans =
% 
%        0     0     0     0     0     0     1     1     1     2
%        0     0     0     1     1     2     0     0     1     0
%        0     1     2     0     1     0     0     1     0     0

arguments
    dim
    n
    type = '='
end

switch type
    case '='
        idx = [];
        if ( dim <= 0 || n < 0 ), return, end
        if ( dim == 1 ), idx = n; return, end
        for k = 0:n-1
            idx1 = tuples(dim-1, n-k, '=');
            idx = [idx [k*ones(1,size(idx1,2)) ; idx1]]; %#ok<AGROW>
        end
        idx = [idx [n ; zeros(dim-1, 1)]];
    case '<='
        idx = tuples(dim+1, n, '=');
        idx = idx(1:dim,:);
    case '<'
        idx = tuples(dim, n-1, '<=');
    otherwise
        error('Type must be either ''='' or ''<=''.');
end

end
