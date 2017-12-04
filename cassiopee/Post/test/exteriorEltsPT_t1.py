# - exteriorElts (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

a = G.cartTetra((0,0,0), (1,1.,1), (20,20,20))
a = C.addVars(a,'Density'); a = C.addVars(a,'centers:cellN')
b = P.exteriorFaces(a)
test.testT(b, 1)

