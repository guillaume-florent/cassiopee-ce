# - triangulateExteriorFaces (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import KCore.test as test

t = C.convertFile2PyTree('boolNG_M1.tp')
t = C.initVars(t, 'F={CoordinateX}')
t = C.initVars(t, 'centers:G={centers:CoordinateX}')
t = C.convertArray2NGon(t)

t = XOR.triangulateExteriorFaces(t)
C.convertPyTree2File(t, 'out.cgns')
test.testT(t, 1)
