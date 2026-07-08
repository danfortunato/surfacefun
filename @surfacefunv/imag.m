function fi = imag(f)
%Imag   Imaginary part of a SURFACEFUNV

fi = f;
fi.components{1} = imag(f.components{1});
fi.components{2} = imag(f.components{2});
fi.components{3} = imag(f.components{3});

end

