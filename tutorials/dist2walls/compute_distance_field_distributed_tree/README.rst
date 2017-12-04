 This tutorial demonstrates how to compute the distance to walls for a distributed pyTree. A distributed pyTree means that zones are distributed on different processors.

          Compute the distance field for a distributed pyTree, by C. Benoit.

The first script is only used to setup the case. Second script performs the computation of wall distance. It is executed on each processor and must be launched with mpirun -np 10 python dist2Walls.py. It loads walls.cgns on every processors, then it loads only the skeleton of the pyTree to be computed. Then a distribution is performed to determine which zone must be loaded on each processor. Since wall distance computation doesnt require any communications once the walls are defined, useCom=0 and fast algorithm are chosen. Then, only certain zones are loaded on processor acording to the distribution using C.readZones. Skeleton tree is then converted to a partial tree, that is a tree where only loaded zones are defined. Dist2Walls are computed on local zones. Finaly, each processor write its zones to a common file using Cmpi.convertPyTree2File.

[Download case script].
[Download python script].
