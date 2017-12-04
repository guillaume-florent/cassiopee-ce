This tutorial demonstrates how to blank a distributed tree, that is a pyTree with zones distributed on different processors.

          Blanking a distributed pyTree, by C. Benoit.

The first script is only used to setup the case. Second script performs blanking. It is executed on each processor and must be launched with mpirun -np 10 python blank.py. It loads walls.cgns on every processors. This file contains the definition of bodies. Then it loads only the skeleton of the pyTree to be blanked. A distribution is performed to determine which zone must be loaded on each processor. Since blanking doesnt require any communications once the walls are defined, useCom=0 and fast algorithm are chosen. Then, only certain zones are loaded on processor acording to the distribution using C.readZones. Skeleton tree is then converted to a partial tree, that is a tree where only loaded zones are defined. Blanking is performed on local zones. Finaly, each processor write its zones to a common file using Cmpi.convertPyTree2File.

[Download case script].
[Download python script].