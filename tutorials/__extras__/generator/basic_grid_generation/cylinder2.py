#!/usr/bin/python
# coding: utf-8

r"""cylinder2 (array)

Create an irregular cylindrical grid (or a portion of cylinder between tetas
and tetae) with ni x nj x nk points, of center-bottom point (xo,yo,zo),
of inner radius R1, outer radius R2, height H and with distributions in r, 
teta, z.
Distributions are arrays defining 1D meshes (x and i varying) 
giving a distribution in [0,1]. Their number of points gives ni, nj, nk.

Parameters:
    (xo,yo,zo) (3-tuple of floats) – coordinates of the starting point
    R1 (float) – value of inner radius
    R2 (float) – value of outer radius
    tetas (float) – start angle (in degree)
    tetae (float) – end angle (in degree)
    H (float) – value of cylinder height
    arrayR (array) – distribution along radius
    arrayTeta (array) – distribution along azimuth
    arrayZ (array) – distribution along height

Returns: a 3D structured mesh
Return type: array or pyTree zone

"""
import Converter as C
import Generator as G
import CPlot

r = G.cart((0., 0., 0.), (0.1, 1., 1.), (11, 1, 1))
teta = G.cart((0., 0., 0.), (0.1, 1., 1.), (11, 1, 1))
z = G.cart((0., 0., 0.), (0.1, 1., 1.), (11, 1, 1))
cyl = G.cylinder2( (0., 0., 0.), 0.5, 1., 360., 0., 10., r, teta, z)
C.convertArrays2File([cyl], "out.plt")
CPlot.display(cyl)


