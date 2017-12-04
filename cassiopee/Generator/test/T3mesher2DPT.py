# - T3mesher2D (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

a = D.circle( (0,0,0), 1, N=50 )
a = C.convertArray2Tetra(a); a = G.close(a)
b = G.T3mesher2D(a, triangulateOnly=0)
t = C.newPyTree(['Base',2]); t[2][1][2].append(b)
C.convertPyTree2File(t, 'out.cgns')
