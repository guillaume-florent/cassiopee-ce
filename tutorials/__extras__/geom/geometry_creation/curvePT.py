#!/usr/bin/env python
# coding: utf-8

r"""curve (pyTree)"""

import Converter.PyTree as C
import Geom.PyTree as D


# User definition of parametric curve
def f(t):
    r"""Parametric function"""
    x = t
    y = t * t + 1
    z = 0.
    return x, y, z

a = D.curve(f)
C.convertPyTree2File(a, 'out.cgns')
