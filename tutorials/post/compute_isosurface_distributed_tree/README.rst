This tutorial demonstrates how to compute an isosurface for a distributed pyTree. A distributed pyTree means that zones are distributed on different processors.

          Compute an isosurface for a distributed pyTree, by C. Benoit.

The first script is only used to setup the case. Second script performs the computation of isosurface. It is executed on each processor and must be launched with mpirun -np 10 python iso.py. The case skeleton is loaded on every processors. Then a distribution is performed to determine which zone must be loaded on each processor. Since isosurface computation doesnt require any communications, useCom=0 and fast algorithm are chosen. Then, only certain zones are loaded on processor acording to the distribution using C.readZones. Skeleton tree is then converted to a partial tree, that is a tree where only loaded zones are defined. P.IsosurfMC is then computed on local zones. Generated isosurfaces are not affected to any processors and must be in order Cmpi.convertPyTree2File to work properly. This is done using Cmpi.setProc.

[Download case script].
[Download python script].