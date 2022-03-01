function connectivity = buildConnectivity(dom)

[nv, nu] = size(dom.x{1});
%corners = [1 1; nv 1; nv nu; 1 nu];
corners = [1 1; nv 1; 1 nu; nv nu];
cidx = sub2ind([nv nu], corners(:,1), corners(:,2));
cidx = cidx(:);

% Node j has (x,y,z) coordinate nodes(j,:)
nodes = zeros(4*length(dom), 3);
j = 0;
for k = 1:length(dom)
    nodes(j+(1:4),:) = [dom.x{k}(cidx) dom.y{k}(cidx) dom.z{k}(cidx)];
    j = j+4;
end
[nodes, ~, elem2node] = unique(nodes, 'rows', 'stable');

% Element i has corner nodes elem2node(i,:)
elem2node = reshape(elem2node, 4, length(dom)).';

% Node j touches elements node2elem{j}
node2elem = cell(length(nodes), 1);
for k = 1:length(dom)
    for j = elem2node(k,:)
        node2elem{j} = [node2elem{j} k];
    end
end

% Edge j has nodes edges(j,:)
edges = zeros(4*length(dom), 2);
j = 0;
for k = 1:length(dom)
    edges(j+1,:) = elem2node(k, [1 2]);
    edges(j+2,:) = elem2node(k, [2 4]);
    edges(j+3,:) = elem2node(k, [4 3]);
    edges(j+4,:) = elem2node(k, [3 1]);
%     edges(j+1,:) = elem2node(k, [1 2]);
%     edges(j+2,:) = elem2node(k, [2 3]);
%     edges(j+3,:) = elem2node(k, [3 4]);
%     edges(j+4,:) = elem2node(k, [4 1]);
    j = j+4;
end
edges = sort(edges, 2);
[edges, ~, elem2edge] = unique(edges, 'rows', 'stable');

% Element i has edges elem2edge(i,:)
elem2edge = reshape(elem2edge, 4, length(dom)).';

% Edge j touches elements edge2elem{j}
edge2elem = cell(length(edges), 1);
for k = 1:length(dom)
    for j = elem2edge(k,:)
        edge2elem{j} = [edge2elem{j} k];
    end
end

for j = 1:length(edge2elem)
    if ( length(edge2elem{j}) < 2 )
        edge2elem{j} = [edge2elem{j} -1];
    end
end

% Element i has neighboring elements elem2elem(i,:)
elem2elem = zeros(length(dom), 4);
elem2elem_all = cell2mat(edge2elem(elem2edge));
for k = 1:length(dom)
    idx = elem2elem_all(k,:) ~= k;
    if ( any(idx) )
        elem2elem(k,:) = elem2elem_all(k, idx);
    else
        elem2elem(k,:) = -1;
    end
end

connectivity = struct();
connectivity.nodes     = nodes;
connectivity.elem2node = elem2node;
connectivity.node2elem = node2elem;
connectivity.edges     = edges;
connectivity.elem2edge = elem2edge;
connectivity.edge2elem = edge2elem;
connectivity.elem2elem = elem2elem;

end
