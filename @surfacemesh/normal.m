function normal = normal(dom)

nx = cell(length(dom), 1);
ny = cell(length(dom), 1);
nz = cell(length(dom), 1);
for k = 1:length(dom)
    nx{k} = dom.facenormals{k}(:,:,1);
    ny{k} = dom.facenormals{k}(:,:,2);
    nz{k} = dom.facenormals{k}(:,:,3);
end

normal = surfacefunv(nx, ny, nz, dom);

end
