# - -
import Converter as C
import Intersector as XOR
import KCore.test as test

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

m=XOR.extractOuterLayers(M1)

test.testA(m,1)
