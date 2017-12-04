This tutorial demonstrates how to perform Chimera assembly of a mesh composed by two overlapping cylinders and a Cartesian background grid.

              Chimera assembly by S. Péron.

First, a cellN variable is created. cellN will be equal to 1 for computed cells, 0 for blanked cells and 2 for interpolated cells.
To blank cells inside bodies defined by the two cylinders, we use the function X.blankCells of the Connector module.
To perform an overlap optimization, we use first the function X.applyBCOverlaps to set the cellN to 2 to cells near overlap borders. Then, functions X.optimizeOverlap and X.maximizeBlankedCells are applied to optimize the overlapping between grids.

[Download python script].
