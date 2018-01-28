#!/usr/bin/env python
# coding: utf-8

r"""point (array)"""

import os
import sys

PY3 = False if sys.version_info[0] < 3 else True

import Geom as D
import Converter as C

a = D.point((0,0,0))
C.convertArrays2File([a], "out.plt")

if PY3:
    os.chmod("out.plt", 0o666)
else:
    os.chmod("out.plt", 0666)