function bb = boundingbox(dom)

minall = @(x) min(cellfun(@(y) min(y,[],'all'), x));
maxall = @(x) max(cellfun(@(y) max(y,[],'all'), x));
bb = [minall(dom.x) maxall(dom.x) minall(dom.y) maxall(dom.y) minall(dom.z) maxall(dom.z)];

end
