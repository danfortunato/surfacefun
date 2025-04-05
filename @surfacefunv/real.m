function fr = real(f)
%REAL   Real part of a SURFACEFUNV

fr = f;
fr.components{1} = real(f.components{1});
fr.components{2} = real(f.components{2});
fr.components{3} = real(f.components{3});

end

