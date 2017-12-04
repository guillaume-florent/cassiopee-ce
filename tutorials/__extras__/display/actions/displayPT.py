# - display (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import CPlot.PyTree
import Transform.PyTree as T

a = G.cart((0,0,0),(1,1,1),(18,28,3))
t = C.newPyTree(['Base',a])

for i in xrange(360):
    t = T.rotate(t, (9, 14, 3.5), (0,0,1), 1.)
    CPlot.PyTree.display(t)
