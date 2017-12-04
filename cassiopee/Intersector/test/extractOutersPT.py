# - extractOuterLayers (pyTree) -
import Converter.PyTree as C
import Intersector.PyTree as XOR

t = C.convertFile2PyTree('boolNG_M1.tp')
t = XOR.extractOuterLayers(t)

C.convertPyTree2File(t, "out.cgns")
