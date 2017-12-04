This tutorial demonstrates how to cut an octree mesh below (or above) X/Y/Z planes.

              Octree cut below two X and Z planes by T. Renaud.

By default, an octree has the same size in all directions (dfar parameter). To cut an octree with respect to X/Y/Z planes, we use the functions T.symetrize and P.selectCells.
T.symetrize is used to ensure that an octree matching boundary is located on the X/Y/Z plane. Then, P.selectCells is used to choose the cells from one or other side of the X/Y/Z plane.

[Download python script].