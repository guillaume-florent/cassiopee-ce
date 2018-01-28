#!/usr/bin/env python
# coding: utf-8

r"""line (pyTree)"""

import Geom.PyTree as D
import Converter.PyTree as C

a = D.line((0, 0, 0), (1, 0, 0))
C.convertPyTree2File(a, 'out.cgns')
