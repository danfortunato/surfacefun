function dom = flipz(dom)

zz = reshape([dom.z{:}], [], 1);
shift = min(zz) + max(zz);

for k = 1:length(dom)
    dom.z{k} = shift - dom.z{k};
end

end
