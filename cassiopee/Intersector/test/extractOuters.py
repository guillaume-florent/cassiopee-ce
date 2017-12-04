# - extractOuterLayers (array) -
import Converter as C
import Intersector as XOR

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

m = XOR.extractOuterLayers(M1)

C.convertArrays2File(m, "out.plt")
