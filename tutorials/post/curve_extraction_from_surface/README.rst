This tutorial demonstrates how to extract a curve from a surface.

              Extraction of curves from a surface sphere by C. Benoit.

To extract a curve defined by x=c on a given surface, we use the function P.isoSurfMC. If a solution is defined on the input surface, it is also defined in the resulting curves. If you want to extract the curves defined by a function of x, y, z, say for instance, x+y = c, you must first defined an intermediate variable v = x+y in the surface by: b = C.initVars(a, {v} = {CoordinateX}+{CoordinateY}) and use v in the extraction.

[Dowload python script].