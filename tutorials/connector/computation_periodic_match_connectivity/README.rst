This tutorial demonstrates how to compute the periodic match connectivity in the case of an azimuthal periodicity.

              Computes the periodic match connectivity of a quarter of cylinder by S. Péron.

Suppose you have a periodic mesh in azymuth and want to compute the periodic match connectivity.
The idea is to duplicate the meshes by rotation and then perform a X.connectMatch on the full mesh. Since zone names in the duplicated mesh are the same as in original mesh, this will result in a correct connectivity. Then, we perform the same thing for the other direction.

[Download python script].
