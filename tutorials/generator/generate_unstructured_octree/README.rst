This tutorial demonstrates how to generate an unstructured octree around a given surface. The body surface is digging a hole into the octree. The created cavity is then filled with tetras. An additional boundary prismatic layer can eventually be added. The resulting mesh is provided as a single polyhedral conforming zone, with boundary conditions set.

              Unstructured octrees (without and with boundary layer) around a cat body.

The input of the script must be a triangulated surface mesh, closed and watertight, and already correctly meshed for CFD.

If hWall is positive, a prismatic boundary layer is generated from the input surface using G.addNormalLayers. By default, the boundary layer consists in two extrusions with different geometric series.

An unstructured octree is generated around the surface (the body surface or the external face of the boundary layer) using G.octree. By default, the snear parameter is automatically calculated with respect to the averaged size of the surface elements. By default, the dfar parameter is equal to 10 times the maximum size of the body. If needed, these parameters can be modified in the script.

If the input case contains a 'Refine' base, then the surfaces defined in this base are used to refine locally the octree using G.adaptOctree.

An offset surface, located to a certain distance of the input surface (the body surface or the external face of the boundary layer) is then computed. This surface diggs a hole in the octree using X.blankCells.

Then the space is filled with tetras using G.tetraMesher.