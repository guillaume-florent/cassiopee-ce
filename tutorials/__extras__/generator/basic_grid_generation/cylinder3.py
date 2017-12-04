#!/usr/bin/python
# coding: utf-8

r"""cylinder3 (array)

Create an irregular cylindrical grid (or a portion of cylinder between 
tetas and tetae) from a xz plane mesh defined by a and a teta distribution
defined by arrayTeta.
Parameters:
    a ([array, list of arrays] or [zone, list of zones, base, pyTree]) – 
       definition of the xz plane mesh
    tetas (float) – start angle (in degree)
    tetae (float) – end angle (in degree)
    arrayTeta (array) – distribution along azimuth

Returns:a 3D structured mesh
Return type:array or pyTree zone

"""
import Generator as G
import Converter as C
import CPlot

teta = G.cart((0., 0., 0.), (0.1, 1., 1.), (11, 1, 1))
xz = G.cart((0.1, 0., 0.), (0.1, 1., 0.2), (20, 1, 30))
cyl = G.cylinder3(xz, 0., 90., teta)
C.convertArrays2File([cyl], 'out.plt')
CPlot.display(cyl)
