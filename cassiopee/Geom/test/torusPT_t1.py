# - torus (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

# trimmed torus
a = D.torus((0,0,0), 5., 2., 0., 120., 0., 90., 100, 36)
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
test.testT(t, 1)

test.writeCoverage(100)
