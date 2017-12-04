This tutorial demonstrates how to compute the distance to walls in a pyTree and dump the distance fields to files for further reading by elsA solver. If the 'cellN' variable (located at centers, cellN=0 for blanked cells, cellN=2 for interpolated cells, cellN=1 for computed cells) is defined in the pyTree, then it is taken into account during the distance field computation. If a cell near the wall is blanked, then the corresponding wall face is not used for the distance computation.

          Compute the distance field for a multiblock structured mesh around a sphere, by S. Péron.

In the script, wall boundaries are first extracted as 2D zones using C.extractBCOfType of Converter module. To compute the wall distance, one can use the walls in nodes, in centers or in extended centers using C.node2Center or C.node2ExtCenter. Finally, fonction DTW.distance2Walls is used to compute the distance located at centers in the multiblock mesh defined as a pyTree, using "ortho" or "mininterf" algorithm. To obtain results similar to the distance computed by elsA Kernel, one must choose walls in centers and "mininterf" algorithm.
The result is the creation of a 'TurbulentDistance' node in the pyTree.
Finally, these nodes can be dumped to files which can be read by elsA Kernel.

[Download python script].