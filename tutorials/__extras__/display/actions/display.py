# - display (array) -
import Generator as G
import CPlot
import Transform as T
import Converter as C

a = G.cart((0,0,0),(1,1,1),(18,28,3))
CPlot.display(a, displayBB=0, mode='mesh')

for i in xrange(360):
    a = T.rotate(a, (9, 14, 3.5), (0,0,1), 1.)
    CPlot.display(a)
