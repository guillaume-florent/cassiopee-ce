#!/usr/bin/env python
# coding: utf-8

r"""Offscreen rendering

This tutorial demonstrates how to use offscreen rendering and
dump an image or a mpeg movie.

Rotating cube by C. Benoit.

CPlot can render in an image file or in a mpeg movie file.
If you are using CPlot on a cluster (with no graphic card),
it must be installed with mesa.
Then, when using CPlot.display, offscreen is triggered with offscreen=1.
If you are using CPlot on a computer with a graphic card,
you must use offscreen=2.
In this tutorial, we consider a Cartesian grid, make it rotate with T.rotate,
and dump image to a file at each step. The movie must be finalized at the end.
A wait time is added before finalizing to be sure that render is finished.

"""

import os

import Generator as G
import Converter as C
import CPlot
import time
import Transform as T

a = G.cart((0, 0, 0), (1, 1, 1), (5, 5, 5))

# dump a PNG
CPlot.display([a], mode=0, posCam=(10, 2, 2), posEye=(2, 2, 2),
              export='export.png',
              exportResolution='600x600', offscreen=2)
time.sleep(1)

# dump a MPEG movie
for i in xrange(200):
    a = T.rotate(a, (2, 2, 2), (0, 0, 1), 1)
    time.sleep(0.1)
    CPlot.display([a], mode=1, export='export.mpeg',
                  exportResolution='600x600', offscreen=2)
time.sleep(2)
CPlot.finalizeExport()

os._exit(1)
