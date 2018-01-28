#!/usr/bin/env python
# coding: utf-8

r"""sphere6 (pyTree)"""

import Geom.PyTree as D
import Converter.PyTree as C

A = D.sphere6((0, 0, 0), 1., 20)
C.convertPyTree2File(A, 'out.cgns')
