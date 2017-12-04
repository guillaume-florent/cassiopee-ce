# - display (array) -
# Offscreen using FBO
import Generator as G
import CPlot
import Transform as T
import Converter as C
import Geom as D
import time

a = D.sphere((0,0,0), 1, N=200)

# Multi images
CPlot.display([a], offscreen=2, bgColor=1, mode=0, meshStyle=2,
              solidStyle=1, posCam=(0,6,0), export='one.png')
for i in xrange(5):
    a = T.rotate(a, (0,0,0), (0,0,1), 1.)
    CPlot.display([a], offscreen=2, bgColor=1, mode=0, meshStyle=2,
              solidStyle=1, posCam=(0,6,0), export='one%d.png'%i)
    CPlot.finalizeExport()
import os; os._exit(0)

# Mpeg Movie
for i in xrange(50):
    a = T.rotate(a, (0,0,0), (0,0,1), 1.)
    CPlot.display([a], bgColor=1, mode=0,
                  solidStyle=1, posCam=(0,6,0), export='export.mpeg',
                  exportResolution='680x600', offscreen=2)
    CPlot.finalizeExport()
CPlot.finalizeExport(1)
import os; os._exit(0)
