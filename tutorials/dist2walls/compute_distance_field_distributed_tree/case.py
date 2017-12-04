#!/usr/bin/python
# coding: utf-8

r"""<TBC>"""

import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D

a = G.cart((-2, -2, -2), (0.04, 0.04, 0.04), (100, 100, 100))
a = T.splitSize(a, N=100000)
b = D.sphere((0, 0, 0), 1., N=100)

t = C.newPyTree(['Base'])
t[2][1][2] += a
C.convertPyTree2File(t, 'in.cgns')
t = C.newPyTree(['Base'])
t[2][1][2] += [b]
C.convertPyTree2File(t, 'walls.cgns')
