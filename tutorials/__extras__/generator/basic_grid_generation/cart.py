#!/usr/bin/python
# coding: utf-8

r"""cart (array)

Generator.cart((xo, yo, zo), (hi, hj, hk), (ni, nj, nk))

Create a structured Cartesian mesh with ni x nj x nk points starting 
from point (xo,yo,zo) and of step (hi,hj,hk).
Parameters:	
    (xo,yo,zo) (3-tuple of floats) – coordinates of the starting point
    (hi,hj,hk) (3-tuple of floats) – values of advancing step in the three
                                     directions
    (ni,nj,nk) (3-tuple of integers) – number of points in each direction

Returns: a 1D, 2D or 3D structured mesh
Return type: array or pyTree zone

"""
import Converter as C
import Generator as G
import CPlot

a = G.cart((0., 0., 0.), (0.1, 0.1, 0.2), (10, 10, 10))
print(type(a))
print(len(a))
C.convertArrays2File(a, 'out.plt')
CPlot.display(a, displayBB=0, mode='mesh')
