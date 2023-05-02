function dom = flipx(dom)

xx = reshape([dom.x{:}], [], 1);
shift = min(xx) + max(xx);

for k = 1:length(dom)
    dom.x{k} = shift - dom.x{k};
end

end
