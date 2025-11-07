%% Simple 2:1 example

n = 16;
[xx, yy] = chebpts2(n, n, [0 1 0 1]);
x = {xx-1, 0.5*xx, 0.5*xx+0.5, 0.5*xx,     0.5*xx+0.5}.';
y = {yy,   0.5*yy, 0.5*yy,     0.5*yy+0.5, 0.5*yy+0.5}.';
z = repmat({zeros(n)}, size(x));
split = repmat({zeros(1, 4)}, size(x));
split{1} = [0 1 0 0];
dom = surfacemesh(x, y, z, split);

f = -1;
bc = 0;
pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, f);
u = L.solve(bc);

clf
plot(u), hold on, plot(dom), hold off
shg

%% Tree example
%  Requirements: https://github.com/danfortunato/treefun

rng(0)
n = 16;
f = treefun2(@(x,y) exp(-50*(x.^2+y.^2)), n);

[x, y] = leafpts(f);
g = 0.1*randnfun2;
z = cell(size(x));
for k = 1:length(x)
    z{k} = g(x{k}, y{k});
end

% We only need the first four neighbor types: left, right, down, up
nei = f.leafNeighbors(1:4,:);
split = cell(length(x), 1);
id = leaves(f);
for k = 1:length(id)
    split{k} = cellfun(@(x) length(x)>1, nei(:,id(k))).';
end

dom = surfacemesh(x, y, z, split);

f = -1;
bc = 0;
pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, f);
u = L.solve(bc);

clf
plot(u), hold on, plot(dom), hold off
shg
