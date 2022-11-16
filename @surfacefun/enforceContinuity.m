function u = enforceContinuity(u)
%ENFORCECONTINUITY   Enforce continuity across patches.

pdo = [];
pdo.b = 1;
Id = surfaceop(u.domain, pdo, u);
u = Id.solve(0);

end
