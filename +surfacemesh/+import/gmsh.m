function dom = gmsh(file)

fid = fopen(file);

line = strip(fgetl(fid));
while ( line == "$Comments" )
    jump_to_end_block(fid, 'Comments');
    line = strip(fgetl(fid));
end

% Read mesh format
if ( line ~= "$MeshFormat" )
    error('Invalid Gmsh file.');
end
format = split(strip(fgetl(fid)));
if ( length(format) ~= 3 )
    error('Invalid Gmsh file.');
end
version = format{1};
if ( ~any(format{2} == '01') )
    error('Invalid Gmsh file.');
end
isascii = ( format{2} == '0' );
if ( ~isascii )
    error('Binary Gmsh files are unsupported.');
end
data_size = str2double(format{3});
jump_to_end_block(fid, 'MeshFormat');

while true
    line = skip_blank_lines(fid);
    if ( feof(fid) )
        break
    end

    if ( length(line) < 2 || line(1) ~= '$' )
        error('Invalid Gmsh file');
    end

    section = strip(line(2:end));

    switch ( section )
        case 'Nodes'
            [nodes, node_tags, node_entities] = read_nodes(fid);
        case 'Elements'
            [triangles, quads] = read_elements(fid, node_tags);
        otherwise
            jump_to_end_block(fid, section);
    end
end

fclose(fid);

if ( ~isempty(triangles) && ~isempty(quads) )
    error('Meshes with both triangles and quads are not currently supported.');
end

if ( ~isempty(triangles) )
    npts = length(triangles{1});
    if ( ~all(cellfun(@length, triangles) == npts) )
        error('Meshes with nonuniform polynomial order are not currently supported.');
    end
    n = (sqrt(8*npts+1)-1) / 2;
    ii = gmsh2trianglepts(n);
    V = equi2tri(n);
    x = cell(length(triangles), 1);
    y = cell(length(triangles), 1);
    z = cell(length(triangles), 1);
    for k = 1:length(triangles)
        xyz = nodes(:,triangles{k});
        xyz = xyz(:,ii);
        x{k} = V * xyz(1,:).';
        y{k} = V * xyz(2,:).';
        z{k} = V * xyz(3,:).';
    end
    dom = surfacemesh(x, y, z, 'tri');
elseif ( ~isempty(quads) )
    npts = length(quads{1});
    if ( ~all(cellfun(@length, quads) == npts) )
        error('Meshes with nonuniform polynomial order are not currently supported.');
    end
    n = sqrt(npts);
    ii = gmsh2quadpts(n);
    V = equi2cheb(n);
    x = cell(length(quads), 1);
    y = cell(length(quads), 1);
    z = cell(length(quads), 1);
    for k = 1:length(quads)
        xyz = nodes(:,quads{k});
        xyz = xyz(:,ii);
        x{k} = V * reshape(xyz(1,:), n, n) * V.';
        y{k} = V * reshape(xyz(2,:), n, n) * V.';
        z{k} = V * reshape(xyz(3,:), n, n) * V.';
    end
    dom = surfacemesh(x, y, z, 'quad');
end

end

function [nodes, node_tags, node_entities] = read_nodes(fid)

% numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
line = fscanf(fid, '%d', 4);
[numEntityBlocks, numNodes, minNodeTag, maxNodeTag] = dealm(line);

nodes         = zeros(3, numNodes);
node_tags     = zeros(1, numNodes);
node_entities = zeros(2, numNodes);

idx = 0;
for k = 1:numEntityBlocks
    % entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodes(size_t)
    line = fscanf(fid, '%d', 4);
    [entityDim, entityTag, parametric, numNodesInBlock] = dealm(line);
    if ( parametric ~= 0 )
        error('Parametric nodes are not implemented.');
    end
    ixx = idx + (1:numNodesInBlock);
    node_tags(ixx) = fscanf(fid, '%d', numNodesInBlock);
    nodes(:,ixx) = reshape(fscanf(fid, '%g', 3*numNodesInBlock), 3, numNodesInBlock);
    node_entities(1,ixx) = entityDim;
    node_entities(2,ixx) = entityTag;
    idx = idx + numNodesInBlock;
end

jump_to_end_block(fid, 'Nodes');

end

function [triangles, quads] = read_elements(fid, node_tags)

element_types = dictionary;
element_types(15) = {{'vertex', 1}};
element_types(1)  = {{'line', 2}};
element_types(8)  = {{'line', 3}};
element_types(26) = {{'line', 4}};
element_types(27) = {{'line', 5}};
element_types(28) = {{'line', 6}};
element_types(62) = {{'line', 7}};
element_types(63) = {{'line', 8}};
element_types(64) = {{'line', 9}};
element_types(65) = {{'line', 10}};
element_types(66) = {{'line', 11}};
element_types(2)  = {{'triangle', 3}};
element_types(9)  = {{'triangle', 6}};
element_types(20) = {{'triangle', 9}};
element_types(21) = {{'triangle', 10}};
element_types(22) = {{'triangle', 12}};
element_types(23) = {{'triangle', 15}};
element_types(25) = {{'triangle', 21}};
element_types(42) = {{'triangle', 28}};
element_types(43) = {{'triangle', 36}};
element_types(44) = {{'triangle', 45}};
element_types(45) = {{'triangle', 55}};
element_types(46) = {{'triangle', 66}};
element_types(3)  = {{'quad', 4}};
element_types(10) = {{'quad', 9}};
element_types(16) = {{'quad', 8}};
element_types(36) = {{'quad', 16}};
element_types(37) = {{'quad', 25}};
element_types(38) = {{'quad', 36}};
element_types(39) = {{'quad', 12}};
element_types(41) = {{'quad', 20}};
element_types(47) = {{'quad', 49}};
element_types(48) = {{'quad', 64}};
element_types(49) = {{'quad', 81}};
element_types(50) = {{'quad', 100}};
element_types(51) = {{'quad', 121}};

inv_tags = zeros(max(node_tags), 1);
inv_tags(node_tags) = 1:length(node_tags);

% numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
line = fscanf(fid, '%d', 4);
[numEntityBlocks, numElements, minElementTag, maxElementTag] = dealm(line);

triangles = {};
quads = {};

for k = 1:numEntityBlocks
    % entityDim(int) entityTag(int) elementType(int) numElementsInBlock(size_t)
    line = fscanf(fid, '%d', 4);
    [entityDim, entityTag, elementType, numElementsInBlock] = dealm(line);
    if ( ~isKey(element_types, elementType) )
        error('Unsupported element type.');
    end
    elem = element_types{elementType};
    [elemShape, numNodesPerElem] = elem{:};
    % d = fromfile(f, c_size_t, int(num_ele * (1 + num_nodes_per_ele))).reshape((num_ele, -1))
    for j = 1:numElementsInBlock
        line = fscanf(fid, '%d', numNodesPerElem+1);
        elementTag = line(1);
        nodes = line(2:end);
        if ( elemShape == "triangle" )
            triangles{end+1} = inv_tags(nodes); %#ok<*AGROW>
        elseif ( elemShape == "quad" )
            quads{end+1} = inv_tags(nodes);
        end
    end

end

end

function ii = gmsh2trianglepts(n)

ii = zeros(n);
ox = 1;
oy = n;
idx = 1;
nshell = ceil(n/3);
for shell = 1:nshell
    ii(oy,ox) = idx; idx = idx+1;
    if ( idx > n*(n+1)/2 ), break, end

    h = n-3*shell+2;
    ii(oy-h,ox) = idx; idx = idx+1;
    ii(oy,ox+h) = idx; idx = idx+1;
    if ( idx > n*(n+1)/2 ), break, end

    o = n*(ox-1)+oy;
    hypot = (o+3*shell-1):(n+1):(o+n*h-1);
    ii(oy-(1:h-1),ox)    = idx + (0:h-2); idx = idx + h-1;
    ii(hypot)            = idx + (0:h-2); idx = idx + h-1;
    ii(oy,ox+(h-1:-1:1)) = idx + (0:h-2); idx = idx + h-1;

    ox = ox+1;
    oy = oy-1;
end

ii = flipud(ii);
ii = ii(:);
ii(ii == 0) = [];

end

function ii = gmsh2quadpts(n)

ii = zeros(n);
ox = n;
oy = n;
idx = 1;
nshell = ceil(n/2);
for shell = 1:nshell
    ii(oy,ox) = idx; idx = idx+1;
    if ( idx > n^2 ), break, end

    h = n-2*shell;
    ii(oy,ox-h-1)     = idx; idx = idx+1;
    ii(oy-h-1,ox-h-1) = idx; idx = idx+1;
    ii(oy-h-1,ox)     = idx; idx = idx+1;
    if ( idx > n^2 ), break, end

    ii(oy,ox-(1:h))        = idx + (0:h-1); idx = idx + h;
    ii(oy-(1:h),ox-h-1)    = idx + (0:h-1); idx = idx + h;
    ii(oy-h-1,ox-(h:-1:1)) = idx + (0:h-1); idx = idx + h;
    ii(oy-(h:-1:1),ox)     = idx + (0:h-1); idx = idx + h;

    ox = ox-1;
    oy = oy-1;
end

ii = flipud(ii);
ii = ii(:);

end

function V = equi2tri(n)
    [x, y] = trianglepts(n);
    [xgmsh, ygmsh] = trianglepts(n, type='linspace');
    K     = koornwinder(n-1, x, y);
    Kgmsh = koornwinder(n-1, xgmsh, ygmsh);
    V = K / Kgmsh;
end

function V = equi2cheb(n)
    xgmsh = linspace(-1, 1, n).';
    cheb2equi = bary(xgmsh, eye(n));
    V = inv(cheb2equi);
end

function jump_to_end_block(fid, block)
    found = false;
    while ~feof(fid)
        line = fgetl(fid);
        if strip(line) == sprintf("$End%s", block)
            found = true;
            break
        end
    end
    if ( ~found )
        error('$%s not closed by $End%s.', block);
    end
end

function line = skip_blank_lines(fid)
    while ~feof(fid)
        line = strip(fgetl(fid));
        if ~isempty(line)
            break
        end
    end
end

function varargout = dealm(x)
    n = size(x, 2);
    if ( n == 1 )
        n = size(x, 1);
        x = x(:).';
    end
    if ( nargout > n )
        error('Too many outputs.');
    end
    varargout = cell(1, nargout);
    for k = 1:nargout
        varargout{k} = x(:,k);
    end
end

function y = quadTensorProductOrder(x, n)

y = zeros(n);

y(1) = x(1);

if ( n > 1 )
    y(n) = x(2);
    y(n^2) = x(3);
    y(n*(n-1)+1) = x(4);
end

for i = 1:n-2
    y(i+1) = x(4+i);
end

for i = 1:n-2
    y(n^2-i) = x(4+2*(n-2)+i);
end

% Interior (recursion)
if ( n > 2 )
    % Reorder interior nodes
    interior = x((4+4*(n-2)+1):end);
    interior = quadTensorProductOrder(interior, n-2);

    % Copy into full node list
    for j = 1:n-2
        y(j*n+1) = x(4+3*(n-2)+n-2-j+1);
        for i = 1:n-2
            y(j*n+i+1) = interior((j-1)*(n-2)+(i-1)+1);
        end
        y((j+1)*n) = x(4+(n-2)+j);
    end
end

end
