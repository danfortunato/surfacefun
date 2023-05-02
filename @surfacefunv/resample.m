function u = resample(u, n)

if ( isempty(u) )
    return
end

u.components{1} = resample(u.components{1}, n);
u.components{2} = resample(u.components{2}, n);
u.components{3} = resample(u.components{3}, n);

end
