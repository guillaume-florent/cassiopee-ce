This tutorial demonstrates how to get a full configuration case from a symmetric case with respect to a plane, keeping solution, boundary conditions and connectivity.

          Symmetrize half a cylinder case with respect to y=0 plane by S. Péron.

To get a full case from a symetric case according to a plane x=a (resp. y=a, z=a), we first symmetrize the mesh using the function T.symetrize. The generated blocks are generally ordered in an indirect way. To get this right, we use T.reorder. Finally, the BCSymmetryPlane boundary condition is removed and replaced by match boundary conditions using X.connectMatch.
As a result, a full case with solution and correct boundary conditions is obtained.

[Download python script].