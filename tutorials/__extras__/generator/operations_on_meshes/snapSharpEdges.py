# - snapSharpEdges (array) -
import Generator as G
import Converter as C
import Geom as D

# polylines with sharp angles
s = D.polyline([(0.2, 0, 0), (1, 1, 0), (2.5, 1, 0), (0.2, 0, 0)])
# Regular cartesian grid
h = 0.1
ni = 30
nj = 20
nk = 1
b = G.cart((-0.5, -0.5, 0), (h, h, 1.), (ni, nj, nk))
b = G.snapSharpEdges(b, [s], h)
C.convertArrays2File([b, s], 'out.plt')
