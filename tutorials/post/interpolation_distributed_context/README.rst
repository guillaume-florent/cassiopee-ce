This tutorial demonstrates how to interpolate the field of a donor mesh to a receiver mesh in a distributed context.

          Interpolation with distributed pyTrees, by C. Benoit.

Function P.extractMesh provides interpolation from one donor mesh to a receiver mesh. In order to make interpolation process faster, it is possible to use it in a distributed context on N processors. First the case must be prepared. We suggest to over split the donor mesh in 3*N parts and to split to the receiver mesh in N parts. When this is done, Pmpi.extractMesh will only get the necessary blocks for interpolation on one processor.

[Download case script].
[Download python script].