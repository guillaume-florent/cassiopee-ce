#!/usr/bin/python
# coding: utf-8

r"""cartNGon (array)

Generator.cartNGon((xo, yo, zo), (hi, hj, hk), (ni, nj, nk))

Create a NGON mesh defined from a regular Cartesian mesh. The initial Cartesian
mesh is defined by ni x nj x nk points starting from point (xo,yo,zo) 
and of step (hi,hj,hk). Type of elements is ‘NGON’.
Parameters:
    (xo,yo,zo) (3-tuple of floats) – coordinates of the starting point
    (hi,hj,hk) (3-tuple of floats) – values of advancing step in the three
                                     directions
    (ni,nj,nk) (3-tuple of integers) – number of points in each direction

Returns:a 1D, 2D or 3D unstructured mesh
Return type:array or pyTree zone

"""
import Generator as G
import Converter as C
import CPlot

a = G.cartNGon((0., 0., 0.), (0.1, 0.1, 0.2), (20, 20, 20))
CPlot.display([a], meshStyle=1)
