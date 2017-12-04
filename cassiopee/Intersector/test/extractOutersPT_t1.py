# - extractOuters (pyTree) -
import Converter.PyTree as C
import Intersector.PyTree as XOR
import KCore.test as test

t = C.convertFile2PyTree('boolNG_M1.tp')
t=XOR.extractOuterLayers(t)

test.testT(t, 1)
