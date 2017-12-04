This tutorial demonstrates how to extract slices from a surface.

              Extraction of slices from a surface sphere by C. Benoit.

In this tutorial, we want to extract slices x=c of width dx. We will make use of boolean operator intersection. We first create a box corresponding to a slice and we compute the boolean intersection of the volume defined by this box and the volume defined by the surface. For this to work properly, the normals of the surface must be externaly oriented.

[Dowload python script].