# - setDoublyDefinedBC (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30)); b[0] = 'fente'
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'kmin',[b],'doubly_defined')
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)

t = C.initVars(t, 'centers:cellN', 1)
t = X.applyBCOverlaps(t)
t = X.setDoublyDefinedBC(t)
C.convertPyTree2File(t, 'out.cgns')
