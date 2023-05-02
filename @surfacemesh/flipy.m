function dom = flipy(dom)

yy = reshape([dom.y{:}], [], 1);
shift = min(yy) + max(yy);

for k = 1:length(dom)
    dom.y{k} = shift - dom.y{k};
end

end
