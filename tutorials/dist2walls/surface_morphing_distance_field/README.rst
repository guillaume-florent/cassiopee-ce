This tutorial demonstrates how to perform morphing between two surfaces using the distance field.

              Cassiopee/Onera morphing by C. Benoit.

To perform morphing between two given surfaces, the idea is to use the distance field. We first compute a uniform Cartesian grid around both given surfaces. We then compute the distance field for first surface (Dist1) and for second surface (Dist2) by using Dist2Walls.distance2Walls on the uniform Cartesian grid. Intermediate distance fields are computed by linear interpolation of the two distance fields Dist1 and Dist2. Intermediate surface is reconstructed by marching cubes (using P.isoSurfMC).

[Dowload case].
[Dowload python script].