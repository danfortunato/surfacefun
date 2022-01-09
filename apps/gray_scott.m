% Gray-Scott equations on a surface

% u_t = e1 lap(u) + b(1-u) - uv^2
% v_t = e2 lap(v) - dv + uv^2

%e1 = 0.00002; e2 = 0.00001;
%b = 0.04; d = 0.1;
e1 = 0.00002; e2 = 0.00001;
%F = 0.036; K = 0.057;
F = 0.04; K = 0.06;
dt = 0.2;

Nu = @(u,v) F*(1-u) - u.*v.^2;
Nv = @(u,v) -(F+K)*v + u.*v.^2;

n = 16;
dom = surfacemesh.sphere(n, 2);
pdo = struct('lap', -dt*e1, 'b', 1);
Lu = surfaceop(dom, pdo);
build(Lu)
pdo = struct('lap', -dt*e2, 'b', 1);
Lv = surfaceop(dom, pdo);
build(Lv)

%%
close all
uinit = surfacefun(@(x,y,z) 1-exp(-80*((x+.05).^2+(y+.02).^2 + (z-1).^2)), dom);
vinit = surfacefun(@(x,y,z)   exp(-80*((x-.05).^2+(y-.02).^2 + (z-1).^2)), dom);

rng(0)
%u = randnfunsphere(0.5); u = u - min2(u); u = u ./ max2(u);
%v = randnfunsphere(0.5); v = v - min2(v); v = v ./ max2(v);
% u = randnfunsphere(0.5); u = u ./ max(max2(u), abs(min2(u)));
% v = randnfunsphere(0.5); v = v ./ max(max2(v), abs(min2(v)));
% uinit = (1-0.01) + 0.01*surfacefun(@(x,y,z) u(x,y,z), dom);
% vinit = 0.01*surfacefun(@(x,y,z) v(x,y,z), dom);
% 
% uinit.vals{1} = 0.5  + 0*uinit.vals{1};
% vinit.vals{1} = 0.25 + 0*vinit.vals{1};
% uinit.vals{2} = 0.5  + 0*uinit.vals{1};
% vinit.vals{2} = 0.25 + 0*vinit.vals{1};
% uinit.vals{3} = 0.5  + 0*uinit.vals{1};
% vinit.vals{3} = 0.25 + 0*vinit.vals{1};
% uinit.vals{4} = 0.5  + 0*uinit.vals{1};
% vinit.vals{4} = 0.25 + 0*vinit.vals{1};

%%
close all
u = uinit;
v = vinit;

plot(v)
%view(0,90)
colorbar
shg
pause

t = 0;
for k = 1:10000
    Lu.rhs = u + dt*Nu(u,v);
    Lv.rhs = v + dt*Nv(u,v);
    u = solve(Lu);
    v = solve(Lv);
    t = t + dt;
    if ( mod(k,100) == 0 )
        k
        plot(v)
        %view(0,90)
        colorbar
        drawnow
        shg
    end
end
