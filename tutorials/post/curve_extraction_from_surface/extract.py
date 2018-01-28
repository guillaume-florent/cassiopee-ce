#!/usr/bin/env python
# coding: utf-8

r"""<TBC>"""

# Extraction de lignes x=cte sur un maillage surfacique
import Converter.PyTree as C
import Post.PyTree as P
import Geom.PyTree as D

# IN: surface+solution
a = D.sphere((0, 0, 0), 1, N=100)
a = C.initVars(a, 'p={CoordinateX}+{CoordinateY}')

# Extraction de lignes x = cste
lines = []
for i in range(10):
        x = i/10.
        p = P.isoSurfMC(a, 'CoordinateX', x)
        lines += p

# Sortie
t = C.newPyTree(['Base'])
t[2][1][2] += lines
C.convertPyTree2File(t, 'out.cgns')
