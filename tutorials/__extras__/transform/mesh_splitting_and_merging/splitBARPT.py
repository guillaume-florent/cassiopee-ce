# - splitBAR (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (50,1,1))
a = C.convertArray2Tetra(a)
a = G.close(a)
B = T.splitBAR(a, 5)
C.convertPyTree2File(B, 'out.cgns')



