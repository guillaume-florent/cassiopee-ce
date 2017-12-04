#!/usr/bin/python
# coding: utf-8

r"""Tutorial 1 for Cassiopee"""

from __future__ import print_function

# import Converter.PyTree as C
import Generator.PyTree as G


a = G.cart(Xo=(0, 0, 0), H=(1, 1, 1), N=(10, 10, 10))
print("Creating a cartesian grid.")
