# - point (pyTree) -

import os
import sys

PY3 = False if sys.version_info[0] < 3 else True

import Geom.PyTree as D
import Converter.PyTree as C

a = D.point((0,0,0))
C.convertPyTree2File([a], "out.cgns")

if PY3:
    os.chmod("out.cgns", 0o666)
else:
    os.chmod("out.cgns", 0666)