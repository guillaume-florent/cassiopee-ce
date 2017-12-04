#!/usr/bin/python
# coding: utf-8

r"""cartTetra (pyTree)

Create an unstructured tetrahedral mesh defined from a Cartesian grid 
of ni x nj x nk points starting from point (xo,yo,zo) and of step (hi,hj,hk).
Type of elements are ‘TRI’ for 2D arrays and ‘TETRA’ for 3D arrays.
Parameters:
    (xo,yo,zo) (3-tuple of floats) – coordinates of the starting point
    (hi,hj,hk) (3-tuple of floats) – values of advancing step in the three directions
    (ni,nj,nk) (3-tuple of integers) – number of points in each direction

Returns:a 1D, 2D or 3D unstructured mesh
Return type:array or pyTree zone

"""
import Generator.PyTree as G
import Converter.PyTree as C
import CPlot.PyTree

a = G.cartTetra((0., 0., 0.), (0.1, 0.1, 0.2), (10, 10, 3))
C.convertPyTree2File(a, 'out.cgns')
CPlot.PyTree.display(a, displayBB=0, mode='mesh', meshStyle=1)