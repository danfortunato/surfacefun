function T = trilattice(n)
%TRILATTICE   Indices for a regular lattice of triangles.
%   T = TRILATTICE(N) constructs an (N-1)^2 x 3 matrix T encoding triangle-to-
%   vertex connectivity of the Delaunay triangulation of a regular lattice of
%   nodes inside the unit triangle. The vertices of the i-th subtriangle have
%   indices given by T(i,:).

% Compute triangle-to-vertex connectivity
ntri = (n-1)^2;
T = zeros(ntri, 3);
k = 1;
colstart = 1;
for i = 0:n-2
    % Add triangles spanning columns i and i+1
    h = n-i-1;
    ss = (colstart+1:colstart+h-1).';
    T(k,:)           = [colstart colstart+1 colstart+1+h];
    T(k+1:k+h-1,:)   = [ss ss+h ss+h+1];
    T(k+h:k+2*h-2,:) = [ss ss+1 ss+h+1];
    k = k+2*h-1;
    colstart = colstart+h+1;
end

end
