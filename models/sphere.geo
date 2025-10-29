SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 1};

Mesh.ElementOrder = 10;
Mesh.HighOrderOptimize = 4;
Mesh.MeshSizeFactor = 1;

Mesh 2;
Save "sphere_tri.msh";
Delete Meshes;

Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 1;
Mesh.MeshSizeFactor = 0.7;
Mesh 2;
Save "sphere_quad.msh";
