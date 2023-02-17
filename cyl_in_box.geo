/****************************************/
// Construct the nodes of the boundary and the interior probe/cylinder
Include "mainbox_nodes.geo";

/****************************************/
// Construct the edges of the boundary and the interior probe/cylinder
Include "mainbox_edges.geo";

/****************************************/
// Construct the triangles or patches
Include "mainbox_tris.geo";
/****************************************/
// Add in tetrahedra representing the coarse elements, which also combines all
// the elements
Include "mainbox_vol.geo";

check = 0;

