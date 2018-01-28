#!/usr/bin/env python
# coding: utf-8

r"""point (array)"""

import Geom as D
import Converter as C

a = D.point((0, 0, 0))
C.convertArrays2File([a], "out.plt")
