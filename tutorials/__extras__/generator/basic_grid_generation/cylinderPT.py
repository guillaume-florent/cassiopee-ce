#!/usr/bin/python
# coding: utf-8

r"""cylinder (pyTree)

Create a regular cylindrical grid (or a portion of cylinder between tetas 
and tetae) with ni x nj x nk points, of center-bottom point (xo,yo,zo), 
of inner radius R1, outer radius R2 and height H. For a direct mesh, 
use tetae < tetas.

Parameters:
    (xo,yo,zo) (3-tuple of floats) – coordinates of the starting point
    R1 (float) – value of inner radius
    R2 (float) – value of outer radius
    tetas (float) – start angle (in degree)
    tetae (float) – end angle (in degree)
    (ni,nj,nk) (3-tuple of integers) – number of points in each direction

Returns:a 3D structured mesh
Return type:array or pyTree zone

"""
import Converter.PyTree as C
import Generator.PyTree as G
import CPlot.PyTree
a = G.cylinder((0., 0., 0.), 0.5, 1., 360., 0., 10., (50, 50, 30))
C.convertPyTree2File(a, 'out.cgns')
CPlot.PyTree.display(a)
