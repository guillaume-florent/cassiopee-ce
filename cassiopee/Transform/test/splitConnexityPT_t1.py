# - splitConnexity (pyTree) -
import KCore.test as test
import Transform.PyTree as T
import Geom.PyTree as D
import Converter.PyTree as C

# Non structure + champs
a = D.text2D("CASSIOPEE")
a = C.addVars(a, 'Density'); a = C.initVars(a, 'centers:cellN', 1)
B = T.splitConnexity(a)
t = C.newPyTree(['Base',2]); t[2][1][2] += B
test.testT(t, 1)
