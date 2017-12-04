# - splitTRI (array) -
import Generator as G
import Converter as C
import Geom as D
import Transform as T

a = D.circle( (0,0,0), 1, N=20 )
a = C.convertArray2Tetra(a)
a = G.close(a)
b = G.T3mesher2D(a)
#C.convertArrays2File([b], 'out.plt')
c = [[9, 25, 27, 30, 29, 28, 34, 38, 0], [29, 23, 19, 20, 24, 29]]
d = T.splitTRI(b, c)
C.convertArrays2File([d[0]], 'out1.plt')
C.convertArrays2File([d[1]], 'out2.plt')
